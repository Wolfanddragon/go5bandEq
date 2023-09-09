package main

import (
	"fmt"
	"io"
	"math"
	"os"

	"go-wav"
	"log"

	"github.com/lxn/walk"
	. "github.com/lxn/walk/declarative"
)

var inPth, outPth string         // in\out path to wav files
var cl, clm, cm, chm, ch float32 // coefs for eq bands

type slide struct {
	Min, Max, Value int
}

// Second order filter
type filter2 struct {
	s2, s1, a1, a2, b0, b1, b2 float64
}

// Filter represents a generic filter
type Filter interface {
	// Next runs filter one step with input u and returns output
	Next(float64) float64
}

func wc(fr float64) float64 {
	return 2 * math.Pi * fr * 1. / 44100
}

func valid(wc float64) bool {
	return .0001 < wc && wc < 3.1415 // reasonable bound
}

func valid2(wc, wd float64) bool {
	return .0001 < wc && wc+.0001 < wd && wd < 3.1415
}

// returns prewarped analog cut-off * dt, given desired cut-off * dt
func prewarp(wc float64) float64 {
	x := wc / 2
	x *= x
	return wc * (15 - x) / (15 - 6*x) // wa * dt
}

func (f *filter2) Next(u float64) float64 {
	t := u - f.a1*f.s1 - f.a2*f.s2 // Direct Form 2
	y := f.b0*t + f.b1*f.s1 + f.b2*f.s2
	f.s2, f.s1 = f.s1, t // shift memory
	return y
}

// NewLowPass2 creates second order Low-Pass filter
func NewLowPass2(wc float64) Filter {
	if !valid(wc) {
		return nil
	}
	var f filter2
	wat := prewarp(wc) // wa * dt
	wat2 := wat * wat
	wat *= 2.82842712474619 // sqrt(8)

	a0 := wat2 + 4 + wat
	f.a1 = 2 * (wat2 - 4) / a0
	f.a2 = (wat2 + 4 - wat) / a0
	f.b0 = wat2 / a0
	f.b2 = f.b0
	f.b1 = 2 * f.b0
	f.s2, f.s1 = 0, 0
	return &f
}

// NewHighPass2 creates second order High-Pass filter
func NewHighPass2(wc float64) Filter {
	if !valid(wc) {
		return nil
	}
	var f filter2
	wat := prewarp(wc) // wa * dt
	wat2 := wat * wat
	wat *= 2.82842712474619 // sqrt(8)

	a0 := wat2 + 4 + wat
	f.a1 = 2 * (wat2 - 4) / a0
	f.a2 = (wat2 + 4 - wat) / a0
	f.b0 = 4 / a0
	f.b2 = f.b0
	f.b1 = -2 * f.b0
	f.s2, f.s1 = 0, 0
	return &f
}

// NewBandPass2 creates second order Band-Pass filter
func NewBandPass2(wc, wd float64) Filter {
	if !valid2(wc, wd) {
		return nil
	}
	var f filter2
	wat := prewarp(wc) // wa * dt
	wbt := prewarp(wd) // wb * dt
	wbwa := wbt * wat  // wb * wa * dt^2
	wd = 2 * (wbt - wat)

	a0 := wbwa + 4 + wd
	f.a1 = 2 * (wbwa - 4) / a0
	f.a2 = (wbwa + 4 - wd) / a0
	f.b0 = wd / a0
	f.b2 = -f.b0
	f.b1, f.s2, f.s1 = 0, 0, 0
	return &f
}

func swapValue(data slide) int {
	res := data.Max - data.Min - data.Value
	log.Println(res)
	return res
}

func performEQ() {
	log.Println("Eh brother, welcome to hell!\n")

	//var low, bnd1, bnd2, bnd3, high float32//coefs for frequencies

	file, _ := os.Open(inPth)
	reader := wav.NewReader(file)
	duration, err := reader.Duration()
	if err != nil {
		fmt.Println("Fucked up duration")
	}

	defer file.Close()

	//	wc = 2 * pi * (desired cutoff in Hz) / (sample rate in Hz) =
	//		(desired cutoff in rad/sec) * (sample period in sec)

	flt1l := NewLowPass2(wc(250)) //creating filters for five eq "sliders"
	flt1r := NewLowPass2(wc(250))
	flt2l := NewBandPass2(wc(250), wc(500))
	flt2r := NewBandPass2(wc(250), wc(500))
	flt3l := NewBandPass2(wc(500), wc(1500))
	flt3r := NewBandPass2(wc(500), wc(1500))
	flt4l := NewBandPass2(wc(1500), wc(3000))
	flt4r := NewBandPass2(wc(1500), wc(3000))
	flt5l := NewHighPass2(wc(5000))
	flt5r := NewHighPass2(wc(5000))

	samplesF := []wav.Sample{}
	for { //here we are gathering the samples to work with
		samples, er := reader.ReadSamples()
		if er == io.EOF { //reading samples until the end of the file
			break
		}

		for _, sample := range samples {
			//fmt.Printf("L/R: %d/%d\n", reader.IntValue(sample, 0), reader.IntValue(sample, 1))
			samplesF = append(samplesF, sample)
		}
	}

	fmt.Println(duration)

	samples1 := []wav.Sample{}        //250
	samples2 := []wav.Sample{}        //250-500
	samples3 := []wav.Sample{}        //500-1500
	samples4 := []wav.Sample{}        //1500-3000
	samples5 := []wav.Sample{}        //5000
	for _, sample := range samplesF { //all samples going through filters
		var sampleN1, sampleN2, sampleN3, sampleN4, sampleN5 wav.Sample
		l := flt1l.Next(float64(reader.IntValue(sample, 0)))
		r := flt1r.Next(float64(reader.IntValue(sample, 1)))
		//fmt.Printf("L/R: %v/%v\n", l, r)

		sampleN1.Values[0] = int(l) //create another array for output
		sampleN1.Values[1] = int(r)
		sampleN2.Values[0] = int(flt2l.Next(float64(reader.IntValue(sample, 0))))
		sampleN2.Values[1] = int(flt2r.Next(float64(reader.IntValue(sample, 1))))
		sampleN3.Values[0] = int(flt3l.Next(float64(reader.IntValue(sample, 0))))
		sampleN3.Values[1] = int(flt3r.Next(float64(reader.IntValue(sample, 1))))
		sampleN4.Values[0] = int(flt4l.Next(float64(reader.IntValue(sample, 0))))
		sampleN4.Values[1] = int(flt4r.Next(float64(reader.IntValue(sample, 1))))
		sampleN5.Values[0] = int(flt5l.Next(float64(reader.IntValue(sample, 0))))
		sampleN5.Values[1] = int(flt5r.Next(float64(reader.IntValue(sample, 1))))

		samples1 = append(samples1, sampleN1)
		samples2 = append(samples2, sampleN2)
		samples3 = append(samples3, sampleN3)
		samples4 = append(samples4, sampleN4)
		samples5 = append(samples5, sampleN5)
	}

	for i := 0; i < len(samples1); i++ { //combining all the bands
		samples1[i].Values[0] = int(cl*float32(samples1[i].Values[0]) +
			clm*float32(samples2[i].Values[0]) +
			cm*float32(samples3[i].Values[0]) +
			chm*float32(samples4[i].Values[0]) + ch*float32(samples5[i].Values[0]))

		samples1[i].Values[1] = int(cl*float32(samples1[i].Values[1]) +
			clm*float32(samples2[i].Values[1]) +
			cm*float32(samples3[i].Values[1]) +
			chm*float32(samples4[i].Values[1]) + ch*float32(samples5[i].Values[1]))
	}

	outfile, e := os.Create(outPth)
	if e != nil {
		fmt.Println("You durak")
	}
	defer outfile.Close()

	var numSamples uint32 = uint32(len(samples1)) //uint32(len(samples))
	var numChannels uint16 = 2
	var sampleRate uint32 = 44100
	var bitsPerSample uint16 = 16

	writer := wav.NewWriter(outfile, numSamples, numChannels, sampleRate, bitsPerSample)

	erer := writer.WriteSamples(samples1)
	if erer != nil {
		fmt.Println("You durak sovsem")
	}
}

func gui() {
	var sll, sllm, slm, slhm, slh *walk.Slider
	var lowEdit, hghEdit, lmEdit, mEdit, hmEdit *walk.NumberEdit
	var iPth, oPth *walk.TextEdit

	dataL := slide{0, 100, 50}
	dataLm := slide{0, 100, 50}
	dataM := slide{0, 100, 50}
	dataHm := slide{0, 100, 50}
	dataH := slide{0, 100, 50}

	MainWindow{
		Title:   "5 band EQ",
		Size:    Size{800, 600},
		MinSize: Size{320, 240},
		Layout:  Grid{},
		Children: []Widget{
			Composite{
				Layout:        Grid{Columns: 5},
				StretchFactor: 4,
				Children: []Widget{
					Label{
						ColumnSpan: 5,
						Text:       "Absolute path to input file:",
					},
					TextEdit{
						AssignTo:      &iPth,
						ColumnSpan:    5,
						CompactHeight: true,
					},
					Label{
						ColumnSpan: 5,
						Text:       "Absolute path to output file:",
					},
					TextEdit{
						AssignTo:      &oPth,
						ColumnSpan:    5,
						CompactHeight: true,
					},
					PushButton{
						ColumnSpan: 5,
						Text:       "Start",
						OnClicked: func() {
							inPth = iPth.Text()
							outPth = oPth.Text()
							cl = float32(dataL.Value) / 100
							clm = float32(dataLm.Value) / 100
							cm = float32(dataM.Value) / 100
							chm = float32(dataHm.Value) / 100
							ch = float32(dataH.Value) / 100
							log.Printf("Low: < %f >\n", cl)
							log.Printf("Low-Mid: < %f >\n", clm)
							log.Printf("Mid: < %f >\n", cm)
							log.Printf("High-Mid: < %f >\n", chm)
							log.Printf("High: < %f >\n", ch)
							log.Printf("%s\n", iPth.Text())
							log.Printf("%s\n", oPth.Text())
							performEQ()
						},
					},
					Label{Text: "0-250 Value"},
					Label{Text: "250-500 Value"},
					Label{Text: "500-1500 Value"},
					Label{Text: "1500-3000 Value"},
					Label{Text: "5000-20000 Value"},

					NumberEdit{
						AssignTo: &lowEdit,
						Value:    float64(swapValue(dataL)),
						OnValueChanged: func() {
							dataL.Value = int(lowEdit.Value())
							dataL.Value = swapValue(dataL)
							sll.SetValue(dataL.Value)
						},
					},
					NumberEdit{
						AssignTo: &lmEdit,
						Value:    float64(swapValue(dataLm)),
						OnValueChanged: func() {
							dataLm.Value = int(lmEdit.Value())
							dataLm.Value = swapValue(dataLm)
							sllm.SetValue(dataLm.Value)
						},
					},
					NumberEdit{
						AssignTo: &mEdit,
						Value:    float64(swapValue(dataM)),
						OnValueChanged: func() {
							dataM.Value = int(mEdit.Value())
							dataM.Value = swapValue(dataM)
							slm.SetValue(dataM.Value)
						},
					},
					NumberEdit{
						AssignTo: &hmEdit,
						Value:    float64(swapValue(dataHm)),
						OnValueChanged: func() {
							dataHm.Value = int(hmEdit.Value())
							dataHm.Value = swapValue(dataHm)
							slhm.SetValue(dataHm.Value)
						},
					},
					NumberEdit{
						AssignTo: &hghEdit,
						Value:    float64(swapValue(dataH)),
						OnValueChanged: func() {
							dataH.Value = int(hghEdit.Value())
							dataH.Value = swapValue(dataH)
							slh.SetValue(dataH.Value)
						},
					},
					Slider{
						RowSpan:        5,
						AssignTo:       &sll,
						MinValue:       dataL.Min,
						MaxValue:       dataL.Max,
						Value:          swapValue(dataL),
						Orientation:    Vertical,
						ToolTipsHidden: true,

						OnValueChanged: func() {
							dataL.Value = sll.Value()
							dataL.Value = swapValue(dataL)
							lowEdit.SetValue(float64(dataL.Value))
						},
					},
					Slider{
						RowSpan:        5,
						AssignTo:       &sllm,
						MinValue:       dataLm.Min,
						MaxValue:       dataLm.Max,
						Value:          swapValue(dataLm),
						Orientation:    Vertical,
						ToolTipsHidden: true,
						OnValueChanged: func() {
							dataLm.Value = sllm.Value()
							dataLm.Value = swapValue(dataLm)
							lmEdit.SetValue(float64(dataLm.Value))
						},
					},
					Slider{
						RowSpan:        5,
						AssignTo:       &slm,
						MinValue:       dataM.Min,
						MaxValue:       dataM.Max,
						Value:          swapValue(dataM),
						Orientation:    Vertical,
						ToolTipsHidden: true,
						OnValueChanged: func() {
							dataM.Value = slm.Value()
							dataM.Value = swapValue(dataM)
							mEdit.SetValue(float64(dataM.Value))
						},
					},
					Slider{
						RowSpan:        5,
						AssignTo:       &slhm,
						MinValue:       dataHm.Min,
						MaxValue:       dataHm.Max,
						Value:          swapValue(dataHm),
						Orientation:    Vertical,
						ToolTipsHidden: true,
						OnValueChanged: func() {
							dataHm.Value = slhm.Value()
							dataHm.Value = swapValue(dataHm)
							hmEdit.SetValue(float64(dataHm.Value))
						},
					},
					Slider{
						RowSpan:        5,
						AssignTo:       &slh,
						MinSize:        Size{Height: 120},
						MinValue:       dataH.Min,
						MaxValue:       dataH.Max,
						Value:          swapValue(dataH),
						Orientation:    Vertical,
						ToolTipsHidden: true,
						OnValueChanged: func() {
							dataH.Value = slh.Value()
							dataH.Value = swapValue(dataH)
							hghEdit.SetValue(float64(dataH.Value))
						},
					},
				},
			},
		},
	}.Run()
}

func main() {
	gui()
}
