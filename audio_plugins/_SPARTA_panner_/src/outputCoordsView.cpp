/*
  ==============================================================================

  This is an automatically generated GUI class created by the Projucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Projucer version: 7.0.7

  ------------------------------------------------------------------------------

  The Projucer is part of the JUCE library.
  Copyright (c) 2020 - Raw Material Software Limited.

  ==============================================================================
*/

//[Headers] You can add your own extra header files here...

//[/Headers]

#include "outputCoordsView.h"


//[MiscUserDefs] You can add your own user definitions and misc code here...
const int sensorEdit_width = 216;
const int sensorEdit_height = 32;
//[/MiscUserDefs]

//==============================================================================
outputCoordsView::outputCoordsView (PluginProcessor* ownerFilter, int _maxNCH, int _currentNCH )
{
    //[Constructor_pre] You can add your own custom stuff here..
    //[/Constructor_pre]

    dummySlider.reset (new juce::Slider ("new slider"));
    addAndMakeVisible (dummySlider.get());
    dummySlider->setRange (0.01, 0.3, 0.001);
    dummySlider->setSliderStyle (juce::Slider::LinearHorizontal);
    dummySlider->setTextBoxStyle (juce::Slider::TextBoxRight, false, 70, 20);
    dummySlider->addListener (this);

    dummySlider->setBounds (-176, 144, 96, 16);


    //[UserPreSize]
    //[/UserPreSize]

    setSize (216, 400);


    //[Constructor] You can add your own custom stuff here..
    setSize (sensorEdit_width, sensorEdit_height*currentNCH);
    hVst = ownerFilter;
    hPan = hVst->getFXHandle();
    maxNCH = _maxNCH ;
    currentNCH =_currentNCH;
    aziSliders =  new std::unique_ptr<Slider>[(unsigned long)maxNCH];
    elevSliders = new std::unique_ptr<Slider>[(unsigned long)maxNCH];
    distSliders =  new std::unique_ptr<Slider>[(unsigned long)maxNCH];

    for( int i=0; i<maxNCH; i++){
        /* create and initialise azimuth sliders */
        distSliders[i].reset(new Slider("new slider"));
        addAndMakeVisible(distSliders[i].get());
        distSliders[i]->setRange(0, 10.0, 0.1);
        distSliders[i]->setValue(panner_getLoudspeakerAzi_deg(hPan, i));
        distSliders[i]->setSliderStyle(Slider::LinearHorizontal);
        distSliders[i]->setTextBoxStyle(Slider::TextBoxRight, false, 70, 20);
        distSliders[i]->setBounds(-25, 8 + i * sensorEdit_height, 86, 16);
        distSliders[i]->addListener(this);

        /* create and initialise azimuth sliders */
        aziSliders[i].reset (new Slider ("new slider"));
        addAndMakeVisible (aziSliders[i].get());
        aziSliders[i]->setRange (0, 360.0, 0.1);
        aziSliders[i]->setValue(panner_getLoudspeakerAzi_deg(hPan, i));
      //  aziSliders[i]->setTextBoxStyle(Slider::TextBoxBelow, true, 70, 20);
        aziSliders[i]->setSliderStyle (Slider::LinearHorizontal);
        aziSliders[i]->setTextBoxStyle (Slider::TextBoxRight, true, 70, 20);
        aziSliders[i]->setBounds(50, 8 + i*sensorEdit_height, 86, 16);
        aziSliders[i]->addListener (this);

        /* create and initialise elevation sliders */
        elevSliders[i].reset (new Slider ("new slider"));
        addAndMakeVisible (elevSliders[i].get());
        elevSliders[i]->setRange (-90.0, 90.0, 0.1);
        elevSliders[i]->setValue(panner_getLoudspeakerElev_deg(hPan, i));
        elevSliders[i]->setSliderStyle (Slider::LinearHorizontal);
        elevSliders[i]->setTextBoxStyle (Slider::TextBoxLeft, false, 70, 20);
        elevSliders[i]->setBounds(155, 8 + i*sensorEdit_height, 86, 16);
        elevSliders[i]->addListener (this);
    }

    sliderHasChanged = true;

	refreshCoords();
	resized();

    //[/Constructor]
}

outputCoordsView::~outputCoordsView()
{
    //[Destructor_pre]. You can add your own custom destruction code here..
    //[/Destructor_pre]

    dummySlider = nullptr;


    //[Destructor]. You can add your own custom destruction code here..
    for( int i=0; i<maxNCH; i++){
        aziSliders[i] = nullptr;
        elevSliders[i] = nullptr;
        distSliders[i] = nullptr;
    }
    delete [] aziSliders;
    delete[] elevSliders;
    delete [] distSliders;
    //[/Destructor]
}

//==============================================================================
void outputCoordsView::paint (juce::Graphics& g)
{
    //[UserPrePaint] Add your own custom painting code here..
    //[/UserPrePaint]

    {
        int x = 88, y = 0, width = 88, height = 2048;
        juce::Colour fillColour1 = juce::Colour (0x21ffffff), fillColour2 = juce::Colour (0x05252a25);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setGradientFill (juce::ColourGradient (fillColour1,
                                             88.0f - 88.0f + x,
                                             128.0f - 0.0f + y,
                                             fillColour2,
                                             176.0f - 88.0f + x,
                                             128.0f - 0.0f + y,
                                             false));
        g.fillRect (x, y, width, height);
    }

    {
        int x = 0, y = 0, width = 88, height = 2048;
        juce::Colour fillColour1 = juce::Colour (0x21ffffff), fillColour2 = juce::Colour (0x05252a25);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setGradientFill (juce::ColourGradient (fillColour1,
                                             88.0f - 0.0f + x,
                                             128.0f - 0.0f + y,
                                             fillColour2,
                                             0.0f - 0.0f + x,
                                             128.0f - 0.0f + y,
                                             false));
        g.fillRect (x, y, width, height);
    }

    {
        int x = 172, y = 0, width = 88, height = 2048;
        juce::Colour fillColour1 = juce::Colour(0x21ffffff), fillColour2 = juce::Colour(0x05252a25);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setGradientFill(juce::ColourGradient(fillColour1,
            88.0f - 0.0f + x,
            128.0f - 0.0f + y,
            fillColour2,
            0.0f - 0.0f + x,
            128.0f - 0.0f + y,
            false));
        g.fillRect(x, y, width, height);
    }

    //[UserPaint] Add your own custom painting code here..
    Colour fillColour = Colours::white;
    g.setColour (fillColour);
    g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Regular"));

    for( int i=0; i<maxNCH; i++){
        /* draw sensor IDs */
        g.setColour (fillColour);
        g.drawText (String(i+1), -30, 5+ i*sensorEdit_height, 33, 23,
                    Justification::centred, true);

        /* draw rectangle around sensor parameter */
        //Colour strokeColour = Colour (0x2370702b);
        //g.setColour (strokeColour);
        g.setColour(Colours::white);
        g.setOpacity(0.15f);
        g.drawRect(0, i * sensorEdit_height, sensorEdit_width/3, sensorEdit_height + 1, 1);
        g.drawRect(sensorEdit_width/3, i * sensorEdit_height, sensorEdit_width/3, sensorEdit_height + 1, 1);
        g.drawRect (sensorEdit_width/3*2, i*sensorEdit_height, sensorEdit_width/3, sensorEdit_height+1, 1);
    }


    //[/UserPaint]
}

void outputCoordsView::resized()
{
    //[UserPreResize] Add your own custom resize code here..
    //[/UserPreResize]

    //[UserResized] Add your own custom resize handling here..
    setSize (sensorEdit_width, sensorEdit_height*currentNCH);
    repaint();
    //[/UserResized]
}

void outputCoordsView::sliderValueChanged (juce::Slider* sliderThatWasMoved)
{
    //[UsersliderValueChanged_Pre]
    //[/UsersliderValueChanged_Pre]

    if (sliderThatWasMoved == dummySlider.get())
    {
        //[UserSliderCode_dummySlider] -- add your slider handling code here..
        //[/UserSliderCode_dummySlider]
    }

    //[UsersliderValueChanged_Post]
    //[/UsersliderValueChanged_Post]
}



//[MiscUserCode] You can add your own definitions of your custom methods or any other code here...

void outputCoordsView::refreshCoords(){
    /* update slider values and limits */
    for( int i=0; i<maxNCH; i++){
        aziSliders[i]->setRange (-360.0, 360.0, 0.1);
        aziSliders[i]->setValue(panner_getLoudspeakerAzi_deg(hPan, i), dontSendNotification);
        elevSliders[i]->setRange (-180.0, 180.0, 0.1);
        elevSliders[i]->setValue(panner_getLoudspeakerElev_deg(hPan, i), dontSendNotification);
        distSliders[i]->setRange(0.0, 10.0, 0.1);
        distSliders[i]->setValue(5, dontSendNotification);
    }
}


//[/MiscUserCode]


//==============================================================================
#if 0
/*  -- Projucer information section --

    This is where the Projucer stores the metadata that describe this GUI layout, so
    make changes in here at your peril!

BEGIN_JUCER_METADATA

<JUCER_COMPONENT documentType="Component" className="outputCoordsView" componentName=""
                 parentClasses="public Component" constructorParams="PluginProcessor* ownerFilter, int _maxNCH, int _currentNCH "
                 variableInitialisers="" snapPixels="8" snapActive="1" snapShown="1"
                 overlayOpacity="0.330" fixedSize="1" initialWidth="176" initialHeight="400">
  <BACKGROUND backgroundColour="465323">
    <RECT pos="88 0 88 2048" fill="linear: 88 128, 176 128, 0=21ffffff, 1=5252a25"
          hasStroke="0"/>
    <RECT pos="0 0 88 2048" fill="linear: 88 128, 0 128, 0=21ffffff, 1=5252a25"
          hasStroke="0"/>
  </BACKGROUND>
  <SLIDER name="new slider" id="4689db34530ab7c7" memberName="dummySlider"
          virtualName="" explicitFocusOrder="0" pos="-176 144 96 16" min="0.01"
          max="0.3" int="0.001" style="LinearHorizontal" textBoxPos="TextBoxRight"
          textBoxEditable="1" textBoxWidth="70" textBoxHeight="20" skewFactor="1.0"
          needsCallback="1"/>
</JUCER_COMPONENT>

END_JUCER_METADATA
*/
#endif


//[EndFile] You can add extra defines here...
//[/EndFile]

