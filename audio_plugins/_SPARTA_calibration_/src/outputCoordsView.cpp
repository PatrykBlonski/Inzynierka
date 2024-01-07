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

    dummyLabel.reset (new juce::Label ("new label"));
    addAndMakeVisible (dummyLabel.get());
   // dummylabel->setRange (0.01, 0.3, 0.001);
  //  dummylabel->setlabelStyle (juce::label::LinearHorizontal);
   // dummylabel->setTextBoxStyle (juce::label::TextBoxRight, false, 70, 20);
  //  dummyLabel->addListener (this);

    dummyLabel->setBounds (-176, 144, 96, 16);


    //[UserPreSize]
    //[/UserPreSize]

    setSize (216, 400);


    //[Constructor] You can add your own custom stuff here..
    setSize (sensorEdit_width, sensorEdit_height*currentNCH);
    hVst = ownerFilter;
    hPan = hVst->getFXHandle();
    maxNCH = _maxNCH ;
    currentNCH =_currentNCH;
    aziLabels =  new std::unique_ptr<Label>[(unsigned long)maxNCH];
    elevLabels = new std::unique_ptr<Label>[(unsigned long)maxNCH];
    distLabels =  new std::unique_ptr<Label>[(unsigned long)maxNCH];

    for( int i=0; i<maxNCH; i++){ //outputCordsView.cpp
        /* create and initialise azimuth labels */
        distLabels[i].reset(new Label("new label"));
        addAndMakeVisible(distLabels[i].get());
     //   distLabels[i]->setRange(0, 10.0, 0.1);
     //   distLabels[i]->setValue(calibration_getLoudspeakerAzi_deg(hPan, i));
      //  distLabels[i]->setlabelStyle(label::LinearHorizontal);
      //  distLabels[i]->setTextBoxStyle(label::TextBoxRight, true, 70, 20);
        distLabels[i]->setText(String(calibration_getLoudspeakerDist_deg(hPan, i)), dontSendNotification);
        distLabels[i]->setBounds(45, 8 + i * sensorEdit_height, 86, 16);
     //   distLabels[i]->addListener(this);

        /* create and initialise azimuth labels */
        aziLabels[i].reset (new Label ("new label"));
        addAndMakeVisible (aziLabels[i].get());
      //  azilabels[i]->setRange (0, 360.0, 0.1);
      //  azilabels[i]->setValue(calibration_getLoudspeakerAzi_deg(hPan, i));
      ////  azilabels[i]->setTextBoxStyle(label::TextBoxBelow, true, 70, 20);
      //  azilabels[i]->setlabelStyle (label::LinearHorizontal);
      //  azilabels[i]->setTextBoxStyle (label::TextBoxRight, true, 70, 20);
        aziLabels[i]->setText(String(calibration_getLoudspeakerAzi_deg(hPan, i)), dontSendNotification);
        aziLabels[i]->setBounds(95, 8 + i*sensorEdit_height, 86, 16);
   //     aziLabels[i]->addListener (this);
        
        /* create and initialise elevation labels */
        elevLabels[i].reset (new Label ("new label"));
        addAndMakeVisible (elevLabels[i].get());
       /* elevlabels[i]->setRange (-90.0, 90.0, 0.1);
        elevlabels[i]->setValue(calibration_getLoudspeakerElev_deg(hPan, i));
        elevlabels[i]->setlabelStyle (label::LinearHorizontal);
        elevlabels[i]->setTextBoxStyle (label::TextBoxLeft, true, 70, 20);*/
        elevLabels[i]->setText(String(calibration_getLoudspeakerElev_deg(hPan, i)), dontSendNotification);
        elevLabels[i]->setBounds(150, 8 + i*sensorEdit_height, 86, 16);

    //    elevLabels[i]->addListener (this);
    }

    labelHasChanged = true;

	refreshCoords();
	resized();

    //[/Constructor]
}

outputCoordsView::~outputCoordsView()
{
    //[Destructor_pre]. You can add your own custom destruction code here..
    //[/Destructor_pre]

  //  dummylabel = nullptr;


    //[Destructor]. You can add your own custom destruction code here..
    for( int i=0; i<maxNCH; i++){
        aziLabels[i] = nullptr;
        elevLabels[i] = nullptr;
        distLabels[i] = nullptr;
    }
    delete [] aziLabels;
    delete[] elevLabels;
    delete [] distLabels;
    //[/Destructor]
}

//==============================================================================
void outputCoordsView::paint (juce::Graphics& g)
{
    //[UserPrePaint] Add your own custom painting code here..
    //[/UserPrePaint]

    {
        int x = 176, y = 0, width = 88, height = 2048;
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

   

    //[UserPaint] Add your own custom painting code here..
    Colour fillColour = Colours::white;
    g.setColour (fillColour);
    g.setFont (Font (15.00f, Font::plain).withTypefaceStyle ("Regular"));

    for( int i=0; i<maxNCH; i++){
        /* draw sensor IDs */
        g.setColour (fillColour);
        g.drawText (String(i+1), -5, 5+ i*sensorEdit_height, 33, 23,
                    Justification::centred, true);

        /* draw rectangle around sensor parameter */
        //Colour strokeColour = Colour (0x2370702b);
        //g.setColour (strokeColour);
        g.setColour(Colours::white);
        g.setOpacity(0.15f);
        g.drawRect(0, i * sensorEdit_height, sensorEdit_width, sensorEdit_height + 1, 1);
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

//void outputCoordsView::labelValueChanged (juce::Label* labelThatWasMoved)
//{
//    //[UserlabelValueChanged_Pre]
//    //[/UserlabelValueChanged_Pre]
//
//    if (labelThatWasMoved == dummyLabel.get())
//    {
//        //[UserlabelCode_dummylabel] -- add your label handling code here..
//        //[/UserlabelCode_dummylabel]
//    }
//
//    //[UserlabelValueChanged_Post]
//    //[/UserlabelValueChanged_Post]
//}



//[MiscUserCode] You can add your own definitions of your custom methods or any other code here...

void outputCoordsView::refreshCoords(){
    /* update label values and limits */
    for( int i=0; i<maxNCH; i++){
       // aziLabels[i]->setRange (-360.0, 360.0, 0.1);
        
        distLabels[i]->setText(String(calibration_getLoudspeakerDist_deg(hPan, i)), dontSendNotification);
        distLabels[i]->setBounds(distLabels[i]->getText().contains("-") ? 40 : 45, 8 + i * sensorEdit_height, 86, 16);

      //  elevLabels[i]->setRange (-180.0, 180.0, 0.1);
        aziLabels[i]->setText(String(calibration_getLoudspeakerAzi_deg(hPan, i)), dontSendNotification);
        aziLabels[i]->setBounds(aziLabels[i]->getText().contains("-") ? 90 : 95, 8 + i * sensorEdit_height, 86, 16);

      //  distLabels[i]->setRange(0.0, 10.0, 0.1);
        elevLabels[i]->setText(String(calibration_getLoudspeakerElev_deg(hPan, i)), dontSendNotification);
        elevLabels[i]->setBounds(elevLabels[i]->getText().contains("-") ? 145 : 150, 8 + i * sensorEdit_height, 86, 16);

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
  <label name="new label" id="4689db34530ab7c7" memberName="dummylabel"
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

