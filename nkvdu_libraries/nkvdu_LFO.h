/*
  ==============================================================================

    nkvdu_LFO.h
    Created: 14 Dec 2018 10:48:27am
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include "nkvdu_memoryless.h"

class simple_lfo
{
public:
    simple_lfo()
    {
        setSampleRate(44100.f);
        _phase = 0.f;
    }
    simple_lfo(float sample_rate)
    {
        setSampleRate(sample_rate);
        _phase = 0.f;
    }
    void setSampleRate(float sample_rate)
    {
        this->sampleRate = sample_rate;
        this->T = 1.f / sample_rate;
    }

    float phasor()
    {
        using namespace nkvdu_memoryless;
        
        _phase += _freq * T;
        _phase = mod_1<float>(_phase);
        return _phase;
    }
    float saw()
    {
        return _phase * 2 - 1;
    }
    float square()
    {
        return _phase > 0.5f ? 1.f : 0.f;
    }
    float tri()
    {
        return _phase > 0.5f ? (_phase - 0.5) * 2 : (1.f - _phase - 0.5) * 2;
    }
    float sine()
    {
        using namespace nkvdu_memoryless;
        float sinewave = shapes.up_cos_LUT(_phase);
        return sinewave;
    }

    // takes in value [0 ... 3] and returns interpolated waveform
    // 0 = tri
    // 1 = saw
    // 2 = square
    // 3 = sine
    float multi()
    {
        using namespace nkvdu_memoryless;
        _wave = clamp<float>(_wave, 0.f, 3.f);
        int selector = (int)(floor(_wave));
        float interpolationAmt = mod_1<float>(_wave);

        switch (selector)
        {
            case 0:
            {
                return linterp<float>(tri(), saw(), interpolationAmt);
                break;
            }
            case 1:
            {
                return linterp<float>(saw(), square(), interpolationAmt);
                break;
            }
            case 2:
            {
                return linterp<float>(square(), sine(), interpolationAmt);
                break;
            }
            case 3:
            {
                return sine();
                break;
            }
        }
        // return 0.f;
    }
    
    float _freq, _wave;
private:
    float sampleRate, T;
    float _phase;
    nkvdu_memoryless::trigTables shapes;
};
