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

#include "calibrationView.h"


//[MiscUserDefs] You can add your own user definitions and misc code here...
const float icon_size = 10.0f;


// LoudspeakerIconComponent.h
#pragma once
#include <JuceHeader.h>



void LoudspeakerIconComponent::setTooltip(const String& newTooltip)
{
    SettableTooltipClient::setTooltip(newTooltip);
}

//[/MiscUserDefs]

//==============================================================================
calibrationView::calibrationView(PluginProcessor* ownerFilter, int _width, int _height)
{
    //[Constructor_pre] You can add your own custom stuff here..
    //[/Constructor_pre]


    //[UserPreSize]
    //[/UserPreSize]

    setSize(460, 460);


    //[Constructor] You can add your own custom stuff here..
    hVst = ownerFilter;
    hPan = hVst->getFXHandle();
    width = _width;
    height = _height;
    //for(int src=0; src<MAX_NUM_INPUTS; src++){
    //    SourceIcons[src].setBounds(width - width*(calibration_getSourceAzi_deg(hPan, src) + 180.0f)/360.f - icon_size/2.0f,
    //                               height - height*(calibration_getSourceElev_deg(hPan, src) + 90.0f)/180.0f - icon_size/2.0f,
    //                               icon_size,
    //                               icon_size);
    //}
    //NSources = calibration_getNumSources(hPan);
    NLoudspeakers = calibration_getNumLoudspeakers(hPan) > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : calibration_getNumLoudspeakers(hPan);
    NLoudspeakersPrev = NLoudspeakers;
    //for(int ls=0; ls<NLoudspeakers; ls++){
    //    LoudspeakerIcons[ls].setBounds(width - width*(calibration_getLoudspeakerAzi_deg(hPan, ls) + 180.0f)/360.f - icon_size/2.0f,
    //                                   height - height*(calibration_getLoudspeakerElev_deg(hPan, ls)+90.0f)/180.0f - icon_size/2.0f,
    //                                   icon_size,
    //                                   icon_size);
    //}

    for (int ls = 0; ls < NLoudspeakers; ls++) {
        LoudspeakerIcons[ls] = new LoudspeakerIconComponent();
        LoudspeakerIcons[ls]->setBounds(width - width * (calibration_getLoudspeakerAzi_deg(hPan, ls) + 180.0f) / 360.f - icon_size / 2.0f,
            height - height * (calibration_getLoudspeakerElev_deg(hPan, ls) + 90.0f) / 180.0f - icon_size / 2.0f,
            icon_size,
            icon_size);


    }
    showInputs = true;
    showOutputs = true;
    sourceIconIsClicked = false;

    //[/Constructor]
}

calibrationView::~calibrationView()
{
    //[Destructor_pre]. You can add your own custom destruction code here..
    //[/Destructor_pre]

    for (int ls = 0; ls < NLoudspeakers; ++ls) {
        delete LoudspeakerIcons[ls];
    }

    //[Destructor]. You can add your own custom destruction code here..
    //[/Destructor]
}

//==============================================================================
void calibrationView::paint(juce::Graphics& g)
{
    //[UserPrePaint] Add your own custom painting code here..
    //[/UserPrePaint]

    {
        int x = 0, y = 0, width = 460, height = 460;
        juce::Colour fillColour1 = juce::Colour(0xff4e4e4e), fillColour2 = juce::Colour(0xff202020);
        juce::Colour strokeColour = juce::Colour(0xff9e9e9e);
        //[UserPaintCustomArguments] Customize the painting arguments here..
        //[/UserPaintCustomArguments]
        g.setGradientFill(juce::ColourGradient(fillColour1,
            248.0f - 0.0f + x,
            0.0f - 0.0f + y,
            fillColour2,
            248.0f - 0.0f + x,
            240.0f - 0.0f + y,
            false));
        g.fillRect(x, y, width, height);
        g.setColour(strokeColour);
        g.drawRect(x, y, width, height, 1);

    }

    //[UserPaint] Add your own custom painting code here..


    /* Draw Grid lines and labels */
    int numGridLinesX = 22;
    int numGridLinesY = 22;
    g.setColour(Colours::white);
    g.setOpacity(0.75f);

    g.drawLine(0.0f, height / 2.0f, width, height / 2.0f, 1.0f);
    g.drawLine(width / 2.0f, 0, width / 2.0f, height, 1.0f);

    for (int i = 0; i <= numGridLinesX; i++) {
        g.setOpacity(0.1f);
        g.drawLine((float)i * width / (float)numGridLinesX, 0, (float)i * width / (float)numGridLinesX, height, 1.0f);
        g.setOpacity(0.75f);
        if (i >= (numGridLinesX - 2) / 2) {
            if ((-20 / 2 + i * (int)20 / (numGridLinesX - 2)) % 5 == 0 || (-20 / 2 + i * (int)20 / (numGridLinesX - 2)) == 0) {
                g.drawText(String((float)((-20 / 2 + i * 20.0f / (numGridLinesX - 2)) / 10.0f)),
                    (float)(i + 2) * width / (float)numGridLinesX - 40, height / 2, 40, 20, Justification::centred, true);
            }
        }
        else if (i < (numGridLinesX - 2) / 2) {
            if ((-20 / 2 + i * (int)20 / (numGridLinesX - 2)) % 5 == 0) {
                g.drawText(String((float)((-20 / 2 + i * 20.0f / (numGridLinesX - 2)) / 10.0f)),
                    (float)i * width / (float)numGridLinesX, height / 2, 40, 20, Justification::centred, true);
            }
        }
    }

    for (int i = 0; i <= numGridLinesY; i++) {
        g.setOpacity(0.1f);
        g.drawLine(0, (float)i * height / (float)numGridLinesY, width, (float)i * height / (float)numGridLinesY, 1.0f);
        g.setOpacity(0.75f);
        if (i < (numGridLinesY - 2) / 2) {
            if ((20 / 2 - i * (int)20 / (numGridLinesY - 2)) % 5 == 0) {
                g.drawText(String((float)((20 / 2 - i * (int)20 / (numGridLinesY - 2)) / 10.0f)),
                    width / 2.0f, (float)i * height / (float)numGridLinesY + 8, 40, 20, Justification::centred, true);
            }
        }
        else if (i > (numGridLinesY - 2) / 2) {
            if ((20 / 2 - i * (int)20 / (numGridLinesY - 2)) % 5 == 0) {
                g.drawText(String((float)((20 / 2 - i * (int)20 / (numGridLinesY - 2)) / 10.0f)),
                    width / 2.0f, (float)(i + 2) * height / (float)numGridLinesY - 20 - 9, 40, 20, Justification::centred, true);
            }
        }
    }
    ///* Draw Grid lines and labels */
    //int numGridLinesX = 22;
    //int numGridLinesY = 22;
    //float stepX = 2.0f / (numGridLinesX - 2); // The step between each grid line for X-axis
    //float stepY = 2.0f / (numGridLinesY - 2); // The step between each grid line for Y-axis
    //g.setColour(Colours::white);
    //g.setOpacity(0.75f);

    //g.drawLine(0.0f, height / 2.0f, width, height / 2.0f, 1.0f);
    //g.drawLine(width / 2.0f, 0, width / 2.0f, height, 1.0f);

    //for (int i = 0; i <= numGridLinesX; i++) {
    //    g.setOpacity(0.1f);
    //    g.drawLine((float)i * width / (float)numGridLinesX, 0, (float)i * width / (float)numGridLinesX, height, 1.0f);
    //    g.setOpacity(0.75f);
    //    float labelValueX = -1.0f + i * stepX;
    //    // Check if the labelValueX is a multiple of 0.1
    //    if (fmod(labelValueX, 0.1f) <= std::numeric_limits<float>::epsilon() || fmod(labelValueX, 0.1f) >= 0.1f - std::numeric_limits<float>::epsilon()) {
    //        g.drawText(String(labelValueX, 1),
    //            (float)i * width / (float)numGridLinesX - 20, height / 2 + 20, 40, 20, Justification::centred, true);
    //    }
    //}

    //for (int i = 0; i <= numGridLinesY; i++) {
    //    g.setOpacity(0.1f);
    //    g.drawLine(0, (float)i * height / (float)numGridLinesY, width, (float)i * height / (float)numGridLinesY, 1.0f);
    //    g.setOpacity(0.75f);
    //    float labelValueY = 1.0f - i * stepY;
    //    // Check if the labelValueY is a multiple of 0.1
    //    if (fmod(labelValueY, 0.1f) <= std::numeric_limits<float>::epsilon() || fmod(labelValueY, 0.1f) >= 0.1f - std::numeric_limits<float>::epsilon()) {
    //        g.drawText(String(labelValueY, 1),
    //            width / 2 - 40, (float)i * height / (float)numGridLinesY - 10, 40, 20, Justification::centred, true);
    //    }
    //}




    if (showOutputs) {
        /* Draw loudspeaker icons */
        for (int ls = 0; ls < NLoudspeakers; ls++) {
            /* icon */
            addAndMakeVisible(LoudspeakerIcons[ls]);
        }
        repaint();
    }
    else {
        for (int ls = 0; ls < NLoudspeakers; ls++) {
            /* icon */
            LoudspeakerIcons[ls]->setVisible(false);
        }
    }

    //if(showInputs){
    //    /* Draw Source icons */
    //    for(int src=0; src<NSources; src++){
    //        /* icon */
    //        //g.setColour(Colour::fromFloatRGBA(1.0-((float)src/(float)NSources), 0.3f, ((float)src/(float)NSources), 1.0f));
    //        g.setColour(Colour::fromFloatRGBA(1.0f, 0.0f, 1.0f, 0.85f));
    //        //setColourGradient(g, (float)src/(float)NSources);
    //        g.setOpacity(0.2f);
    //        g.fillEllipse(SourceIcons[src].expanded(8.0f,8.0f));
    //        g.setOpacity(0.4f);
    //        g.fillEllipse(SourceIcons[src].expanded(4.0f, 4.0f));
    //        g.setOpacity(0.85f);
    //        g.fillEllipse(SourceIcons[src]);
    //        /* icon ID */
    //        g.setColour(Colours::white);
    //        g.setOpacity(0.9f);
    //        g.drawText(String(src+1), SourceIcons[src].expanded(10.0f, 0.0f), Justification::centred, true); // .translated(icon_size, -icon_size)
    //    }
    //}
    //[/UserPaint]
}

void calibrationView::resized()
{
    //[UserPreResize] Add your own custom resize code here..
    //[/UserPreResize]

    //[UserResized] Add your own custom resize handling here..
    //[/UserResized]
}

void calibrationView::mouseDown(const juce::MouseEvent& e)
{
    //[UserCode_mouseDown] -- Add your code here...
    for (int i = 0; i < NSources; i++) {
        Rectangle<int> icon_int;
        icon_int.setBounds(SourceIcons[i].getX(),
            SourceIcons[i].getY(),
            SourceIcons[i].getWidth(),
            SourceIcons[i].getHeight());
        if (icon_int.expanded(4, 4).contains(e.getMouseDownPosition())) {
            sourceIconIsClicked = true;
            indexOfClickedSource = i;
            break;
        }
    }
    //[/UserCode_mouseDown]
}

void calibrationView::mouseDrag(const juce::MouseEvent& e)
{
    //[UserCode_mouseDrag] -- Add your code here...
    if (sourceIconIsClicked) {
        Point<float> point;
        point.setXY((float)e.getPosition().getX() - icon_size / 2.0f, (float)e.getPosition().getY() - icon_size / 2.0f);
        calibration_setSourceAzi_deg(hPan, indexOfClickedSource,
            ((width - (point.getX() + icon_size / 2.0f)) * 360.0f) / width - 180.0f);
        calibration_setSourceElev_deg(hPan, indexOfClickedSource,
            ((height - (point.getY() + icon_size / 2.0f)) * 180.0f) / height - 90.0f);
    }

    //[/UserCode_mouseDrag]
}

void calibrationView::mouseUp(const juce::MouseEvent& e)
{
    //[UserCode_mouseUp] -- Add your code here...
    sourceIconIsClicked = false;
    //[/UserCode_mouseUp]
}



//[MiscUserCode] You can add your own definitions of your custom methods or any other code here...
void calibrationView::refreshPanView()
{
    float x;
    float y;

    //for(int src=0; src<MAX_NUM_INPUTS; src++){
    //    SourceIcons[src].setBounds(width - width*(calibration_getSourceAzi_deg(hPan, src) + 180.0f)/360.f - icon_size/2.0f,
    //                               height - height*(calibration_getSourceElev_deg(hPan, src) + 90.0f)/180.0f - icon_size/2.0f,
    //                               icon_size,
    //                               icon_size);
    //}
    //NSources = calibration_getNumSources(hPan);
    NLoudspeakers = calibration_getNumLoudspeakers(hPan) > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : calibration_getNumLoudspeakers(hPan);

    //TEST

    //calibration_setLoudspeakerAzi_deg(hPan, 0, 66.0f);
    //calibration_setLoudspeakerDist_deg(hPan, 0, 5.0f);
    //calibration_setLoudspeakerElev_deg(hPan, 0, 2.0f);

    //calibration_setLoudspeakerAzi_deg(hPan, 1, 31.0f);
    //calibration_setLoudspeakerDist_deg(hPan, 1, 3.0f);
    //calibration_setLoudspeakerElev_deg(hPan, 1, 3.0f);

    //calibration_setLoudspeakerAzi_deg(hPan, 2, 98.0f);
    //calibration_setLoudspeakerDist_deg(hPan, 2, 9.0f);
    //calibration_setLoudspeakerElev_deg(hPan, 2, 2.0f);
    if (NLoudspeakersPrev != NLoudspeakers) {
        if (NLoudspeakers > NLoudspeakersPrev) {
            for (int ls = NLoudspeakersPrev; ls < NLoudspeakers; ls++) {
                LoudspeakerIcons[ls] = new LoudspeakerIconComponent();
                LoudspeakerIcons[ls]->setBounds(width - width * (calibration_getLoudspeakerAzi_deg(hPan, ls) + 180.0f) / 360.f - icon_size / 2.0f,
                    height - height * (calibration_getLoudspeakerElev_deg(hPan, ls) + 90.0f) / 180.0f - icon_size / 2.0f,
                    icon_size,
                    icon_size);
            }
            NLoudspeakersPrev = NLoudspeakers;
        }
        else {
            for (int ls = NLoudspeakers; ls < NLoudspeakersPrev; ls++) {
                delete LoudspeakerIcons[ls];
            }
            NLoudspeakersPrev = NLoudspeakers;
        }
    }

    for (int ls = 0; ls < NLoudspeakers; ls++) {
        calculateCoordinates(calibration_getLoudspeakerDist_plot(hPan, ls) * 10.0f, calibration_getLoudspeakerAzi_deg(hPan, ls), &x, &y);
        LoudspeakerIcons[ls]->setBounds(width - width * (11 - x) / 22.0f - icon_size / 2.0f,
            height - height * (11 + y) / 22.0f - icon_size / 2.0f,
            icon_size,
            icon_size);
        LoudspeakerIcons[ls]->setTooltip(String(calibration_getLoudspeakerDist_deg(hPan, ls)));
    }
    repaint();
}
//[/MiscUserCode]


//==============================================================================
#if 0
/*  -- Projucer information section --

    This is where the Projucer stores the metadata that describe this GUI layout, so
    make changes in here at your peril!

BEGIN_JUCER_METADATA

<JUCER_COMPONENT documentType="Component" className="calibrationView" componentName=""
                 parentClasses="public Component" constructorParams="PluginProcessor* ownerFilter, int _width, int _height"
                 variableInitialisers="" snapPixels="8" snapActive="1" snapShown="1"
                 overlayOpacity="0.330" fixedSize="1" initialWidth="460" initialHeight="460">
  <METHODS>
    <METHOD name="mouseDown (const MouseEvent&amp; e)"/>
    <METHOD name="mouseDrag (const MouseEvent&amp; e)"/>
    <METHOD name="mouseUp (const MouseEvent&amp; e)"/>
  </METHODS>
  <BACKGROUND backgroundColour="323e44">
    <RECT pos="0 0 480 240" fill="linear: 248 0, 248 240, 0=ff4e4e4e, 1=ff202020"
          hasStroke="1" stroke="1, mitered, butt" strokeColour="solid: ff9e9e9e"/>
  </BACKGROUND>
</JUCER_COMPONENT>

END_JUCER_METADATA
*/
#endif


//[EndFile] You can add extra defines here...
//[/EndFile]

