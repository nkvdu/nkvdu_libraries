/*
  ==============================================================================
    nkvdu_filters.h
    Created: 6 Dec 2018 12:12:44am
    Author:  Nicholas Solem
  ==============================================================================
*/

#pragma once
#include "nkvdu_memoryless.h"

namespace nkvdu_filters {
class filter_abstract
{
public:
    virtual ~filter_abstract() = 0;
    //==============================================================================
    virtual void setSampleRate(float sample_rate)
    {
        this->sampleRate = sample_rate;
        this->T = 1.f / sample_rate;
    }  
    virtual void clear() = 0;
    virtual void updateCutoff(float cutoff_target, float oneOverBlockSize) = 0;
    virtual void updateResonance(float res_target, float oneOverBlockSize) = 0;
    float cutoff_to_g(float cutoff)
    {
        //return cutoff * this->T * 0.5;    // no prewarp
        //return tan(cutoff * T * PI);      // expensive prewarp
        using namespace nkvdu_memoryless;
        
        if ((trig.tan_table) != NULL)
        {
            //float wc = TWOPI * cutoff;
            return (float)this->trig.tan_LUT(cutoff * T / 2.f);
        }
        else
            return 0.f;
    }
    //==============================================================================
    nkvdu_memoryless::trigTables trig;
    float sampleRate, T; // why can't these be static?
    float _oneOverBlockSize, _cutoffTarget, _resonanceTarget;
    float z1;
};
    
inline filter_abstract::~filter_abstract() { }
//==================================================================================


    
class onePole   :   public filter_abstract
{
public:
    onePole()
    {
        this->z1 = 0.0f;
        setSampleRate(44100.f);
    }
    onePole(float sample_rate)
    {
        this->z1 = 0.0f;
        setSampleRate(sample_rate);
    }
    void clear()
    {
        y_n = v_n = z1 = 0.f;
    }
    //==============================================================================
    void updateCutoff()
    {
        if (_cutoffTarget != w_c)
            this->w_c += (_cutoffTarget - this->w_c) * _oneOverBlockSize;
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
        this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
    }

    void updateResonance(float res_target, float oneOverBlockSize)
    {/*no resonance for onepole*/}
    //==============================================================================
    float tpt_lp(float input)
    {
        g = cutoff_to_g(this->w_c);
        v_n = (input - z1) * g / (1 + g);
        y_n = v_n + z1;
        z1 = y_n + v_n;
        // cliptest
        return nkvdu_memoryless::clamp<float>(y_n, -1000.f, 1000.f);
    }
    float tpt_lp(float input, float cutoff)
    {
        g = cutoff_to_g(cutoff);
        v_n = (input - z1) * g / (1 + g);
        y_n = v_n + z1;
        z1 = y_n + v_n;
            // cliptest
        return nkvdu_memoryless::clamp<float>(y_n, -1000.f, 1000.f);
    }
    float tpt_hp(float input)
    {
        return input - tpt_lp(input);
    }
    float tpt_hp(float input, float cutoff)
    {
        return input - tpt_lp(input, cutoff);
    }
    //==============================================================================
private:
    float g, v_n, y_n;
    float w_c;
};
//==============================================================================

// NOTHING SO FAR.
class onePole_nonlinear   :   public onePole
{
public:
/*    onePole_nonlinear()
    {
        this->z1 = 0.0f;
        setSampleRate(44100.f);
    }
    onePole_nonlinear(float sample_rate)
    {
        this->z1 = 0.0f;
        setSampleRate(sample_rate);
    }
 */
private:
};

class fourPole_LP_linear    :   public filter_abstract
{
public:
    /*fourPole_LP_linear()
    {
        setSampleRate(44100.f);
    }
    fourPole_LP_linear(float sample_rate)
    {
        setSampleRate(sample_rate);
    }*/
    fourPole_LP_linear()
        :   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
    {
        setSampleRate(44100.f);
        H1.setSampleRate(44100.f);
        H2.setSampleRate(44100.f);
        H3.setSampleRate(44100.f);
        H4.setSampleRate(44100.f);
    }
    fourPole_LP_linear(float sample_rate)
        :   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
    {       
        setSampleRate(sampleRate);
        H1.setSampleRate(sampleRate);
        H2.setSampleRate(sampleRate);
        H3.setSampleRate(sampleRate);
        H4.setSampleRate(sampleRate);
    }
    void initialize(float sample_rate)
    {
        this->clear();
        setSampleRate(sampleRate);
        H1.setSampleRate(sampleRate);
        H2.setSampleRate(sampleRate);
        H3.setSampleRate(sampleRate);
        H4.setSampleRate(sampleRate);
    }
    void updateOneOverBlockSize(float oneOverBlockSize)
    {
        this->_oneOverBlockSize = oneOverBlockSize;
        H1._oneOverBlockSize = oneOverBlockSize;
        H2._oneOverBlockSize = oneOverBlockSize;
        H3._oneOverBlockSize = oneOverBlockSize;
        H4._oneOverBlockSize = oneOverBlockSize;
    }
    void clear()
    {
        H1.clear();
        H2.clear();
        H3.clear();
        H4.clear();
        s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
    }
    void updateCutoff()
    {
        float local_cutoff = this->w_c;
        if (_cutoffTarget != local_cutoff)
        {
            this->w_c += (_cutoffTarget - local_cutoff) * _oneOverBlockSize;
            H1._cutoffTarget = local_cutoff;
            H2._cutoffTarget = local_cutoff;
            H3._cutoffTarget = local_cutoff;
            H4._cutoffTarget = local_cutoff;
            H1.updateCutoff();
            H2.updateCutoff();
            H3.updateCutoff();
            H4.updateCutoff();
        }
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
        this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
        H1.updateCutoff(cutoff_target, oneOverBlockSize);
        H2.updateCutoff(cutoff_target, oneOverBlockSize);
        H3.updateCutoff(cutoff_target, oneOverBlockSize);
        H4.updateCutoff(cutoff_target, oneOverBlockSize);
    }
    void updateResonance()
    {
        if (_resonanceTarget != k)
            this->k += (_resonanceTarget - this->k) * _oneOverBlockSize;
    }
    void updateResonance(float res_target, float oneOverBlockSize)
    {
        //this->q += (res_target - this->q) * oneOverBlockSize;
        this->k += (res_target - this->k) * oneOverBlockSize;
    }

    float tpt_fourpole(float input)
    {
        g = cutoff_to_g(this->w_c);
        G = g * g * g * g;
        s1 = H1.z1;
        s2 = H2.z1;
        s3 = H3.z1;
        s4 = H4.z1;
        S = g*g*g*s1 + g*g*s2 + g*s3 + s4;

        u_n = ((input) * (1 + k) - k * S) / (1 + k * G);

        y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));

        return y_n;
    }

    float tpt_fourpole(float input, float cutoff)
    {
        g = cutoff_to_g(cutoff);
        G = g * g * g * g;
        s1 = H1.z1;
        s2 = H2.z1;
        s3 = H3.z1;
        s4 = H4.z1;
        S = g*g*g*s1 + g*g*s2 + g*s3 + s4;

        u_n = ((input) * (1 + k) - k * S) / (1 + k * G);

        y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));

        return y_n;
    }
    
    float _resonanceTarget;
private:
    onePole H1, H2, H3, H4;
    float u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
    float w_c, q;
};

// so far quite tame, since the only nonlinearity is in the feedback path. 
// TODO: convert each onepole into a nonlinear onepole
class fourPole_LP_nonlinear    :   public filter_abstract
{
public:
    /*fourPole_LP_linear()
     {
     setSampleRate(44100.f);
     }
     fourPole_LP_linear(float sample_rate)
     {
     setSampleRate(sample_rate);
     }*/
    fourPole_LP_nonlinear()
    :   iters(16),
    u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
    {
        setSampleRate(44100.f);
        H1.setSampleRate(44100.f);
        H2.setSampleRate(44100.f);
        H3.setSampleRate(44100.f);
        H4.setSampleRate(44100.f);
    }
    fourPole_LP_nonlinear(float sample_rate)
    :   iters(16),
    u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
    {
        setSampleRate(sampleRate);
        H1.setSampleRate(sampleRate);
        H2.setSampleRate(sampleRate);
        H3.setSampleRate(sampleRate);
        H4.setSampleRate(sampleRate);
    }
    void initialize(float sample_rate)
    {
        this->clear();
        setSampleRate(sampleRate);
        H1.setSampleRate(sampleRate);
        H2.setSampleRate(sampleRate);
        H3.setSampleRate(sampleRate);
        H4.setSampleRate(sampleRate);
    }
    void updateOneOverBlockSize(float oneOverBlockSize)
    {
        this->_oneOverBlockSize = oneOverBlockSize;
        H1._oneOverBlockSize = oneOverBlockSize;
        H2._oneOverBlockSize = oneOverBlockSize;
        H3._oneOverBlockSize = oneOverBlockSize;
        H4._oneOverBlockSize = oneOverBlockSize;
    }
    void clear()
    {
        H1.clear();
        H2.clear();
        H3.clear();
        H4.clear();
        s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
    }
    void updateCutoff()
    {
        float local_cutoff = this->w_c;
        if (_cutoffTarget != local_cutoff)
        {
            this->w_c += (_cutoffTarget - local_cutoff) * _oneOverBlockSize;
            H1._cutoffTarget = local_cutoff;
            H2._cutoffTarget = local_cutoff;
            H3._cutoffTarget = local_cutoff;
            H4._cutoffTarget = local_cutoff;
            H1.updateCutoff();
            H2.updateCutoff();
            H3.updateCutoff();
            H4.updateCutoff();
        }
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
        this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
        H1.updateCutoff(cutoff_target, oneOverBlockSize);
        H2.updateCutoff(cutoff_target, oneOverBlockSize);
        H3.updateCutoff(cutoff_target, oneOverBlockSize);
        H4.updateCutoff(cutoff_target, oneOverBlockSize);
    }
    void updateResonance()
    {
        if (_resonanceTarget != k)
            this->k += (_resonanceTarget - this->k) * _oneOverBlockSize;
    }
    void updateResonance(float res_target, float oneOverBlockSize)
    {
        //this->q += (res_target - this->q) * oneOverBlockSize;
        this->k += (res_target - this->k) * oneOverBlockSize;
    }
    
    float tpt_fourpole(float input)
    {
        g = cutoff_to_g(this->w_c);
        G = g * g * g * g;
        s1 = H1.z1;
        s2 = H2.z1;
        s3 = H3.z1;
        s4 = H4.z1;
        S = g*g*g*s1 + g*g*s2 + g*s3 + s4;
        
        u_n = y_n;  // initial estimation
        for (int n = 0; n < iters; n++)
            u_n = input - k * (G * tables.tanh_LUT(u_n) + S);
        
        //u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
        
        y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
        
        return y_n;
    }
    
    float tpt_fourpole(float input, float cutoff)
    {
        g = cutoff_to_g(cutoff);
        G = g * g * g * g;
        s1 = H1.z1;
        s2 = H2.z1;
        s3 = H3.z1;
        s4 = H4.z1;
        S = g*g*g*s1 + g*g*s2 + g*s3 + s4;
        
        u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
        
        y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
        
        return y_n;
    }
    
    float _resonanceTarget;
private:
    onePole H1, H2, H3, H4;
    nkvdu_memoryless::trigTables tables;
    int iters;
    float u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
    float w_c, q;
};

class svf_prototype
{
public:
    float lp() {return _outputs.lp;}
    float bp() {return _outputs.bp;}
    float hp() {return _outputs.hp;}
    float np() {return _outputs.np;}

protected:
    struct outputs
    {
        float lp, bp, hp, np;
    } _outputs = { 0.f, 0.f, 0.f, 0.f };
    struct state
    {
        float lp, bp;
    } _state = { 0.f, 0.f };
};

// linear state variable filter using 'naive' integrators (i.e., Euler backward difference integration)
class svf_lin_naive     :   public filter_abstract, svf_prototype
{
public:
    //==============================================================================
    svf_lin_naive() 
    :   w_c(200.f), R(1.f), resonance(1.f)
    {
        this->sampleRate = 44100.f;
        this->T = 1.f / this->sampleRate;
        clear();
    }
    void clear()
    {
        _outputs = {0.f, 0.f, 0.f, 0.f };
        _state = { 0.f, 0.f };
    }
    void setCutoff(float wc)
    {
        this->w_c = wc;
        this->_cutoffTarget = w_c;
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
        this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
    }
    void updateCutoff()
    {
        this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
    }
    void setResonance(float res)
    {
        this->resonance = res;
        this->R = 1.f / res;
    }
    void updateResonance(float res_target, float oneOverBlockSize)
    {
        if (res_target > 0.9f)
            res_target = 0.9f;

        this->resonance += (res_target - this->resonance) * oneOverBlockSize;
        this->R = 1.f / res_target;
    }
    void updateResonance()
    {
        if (this->_resonanceTarget > 0.9f)
            _resonanceTarget = 0.9f;

        this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
        this->R = 1.f / _resonanceTarget;
    }
    /*
    void filter(float input)
    {
        using namespace nkvdu_memoryless;
        _outputs.lp = _state.lp + this->w_c * this->T * _state.bp;
        _outputs.hp = clamp<float>((input - (2 * R * _state.bp) - _state.lp), -10.f, 10.f);
        _outputs.bp = _state.bp + this->w_c * this->T * _outputs.hp;
        
        // update state
        _state.bp = clamp<float>(_outputs.bp, -10.f, 10.f);
        _state.lp = clamp<float>(_outputs.lp, -10.f, 10.f);
    }
    */
    void filter(float input)
    {// Erbe version
        float c, d;
        c = 2.f * sin(PI * w_c * T);
        d = 2.f * (1.f - pow(resonance, 0.25f));

        if (c > 0.5f) c = 0.5f;
        if (d > 2.0f) d = 2.f;
        if (d > (2.f/c - (c * 0.5f)))
            d = 2.f/c - (c * 0.5f);

        _outputs.np = input - (d * _outputs.bp);
        _outputs.lp = _outputs.lp + (c * _outputs.bp);
        _outputs.hp = _outputs.np - _outputs.lp;
        _outputs.bp = _outputs.bp + (c * _outputs.hp);
    }
private:

    float w_c, R, resonance;

/*
    nkvdu_memoryless::trigTables trig;
    float sampleRate, T; // why can't these be static?
    float _oneOverBlockSize, _cutoffTarget, _resonanceTarget;
    float z1;
*/
};

/*
 nonlinear state-variable filter using fourth-order runge-kutta
 y[n+1] = y[n] + 1/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
 h = (1 / fs) / oversample_factor

 k_1 = h*f(t_n, y_n)
 k_2 = h*f(t_n + h/2, y_n + k_1/2)
 k_3 = h*f(t_n + h/2, y_n + k_2/2)
 k_4 = h*f(t_n + h, y_n + k_3)
*/
class svf_nl_rk     :   public filter_abstract, public svf_prototype
{
public:
    svf_nl_rk() 
    : _oversample_factor(1), h(0.000022675736961f),
     w_c(200.f), R(1.f), resonance(1.f)
    {
        this->sampleRate = 44100.f;
        for (int i = 0; i < 2; i++)
        {
            deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
        }
    }

    void setSampleRate(float sample_rate)
    {
        this->sampleRate = sample_rate;
        this->T = 1.f / sample_rate;
        this->h = 1.f / (_oversample_factor * sampleRate);
    }  

    void set_oversample(int oversample_factor)
    {
        _oversample_factor = oversample_factor;
        h = 1.f / (oversample_factor * sampleRate);
    }
    void clear()
    {
        for (int i = 0; i < 2; i++)
        {
            deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
            _outputs = {0.f, 0.f, 0.f, 0.f };
            _state = { 0.f, 0.f };
        }
    }

    void setCutoff(float wc)
    {
        this->w_c = wc;
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
        this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
    }
    void updateCutoff()
    {
        this->w_c += (_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
    }
    void setResonance(float res)
    {
        if (res < 0.5f) res = 0.5f;
        this->resonance = res;
        this->R = 1.f / res;
    }
    void updateResonance(float res_target, float oneOverBlockSize)
    {
        if (res_target < 0.5f) res_target = 0.5f;
        this->resonance += (res_target - this->resonance) * oneOverBlockSize;
        this->R = 1.f / this->resonance;
    }
    void updateResonance()
    {
        if (this->_resonanceTarget < 0.5f) this->_resonanceTarget = 0.5f;
        this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
        this->R = 1.f / this->resonance;
    }
    void filter(float input)
    {
        using namespace nkvdu_memoryless;
        float hp, np;
        // overwritten states. [0] is bp, [1] is lp.
        float tempstate[2];
        
        for (int iter = 0; iter < _oversample_factor; iter++)
        {
            np = input - 2 * R * _state.bp;
            hp = np - _state.lp;
            deriv1[0] = h * TWOPI * w_c * trig.tanh_LUT(hp);
            deriv1[1] = h * TWOPI * w_c * trig.tanh_LUT(_state.lp);

            // I THINK THE PROBLEM IS WITH SCALING, ORDER OF OPERATIONS REGARDING H
            tempstate[0] = _state.bp + deriv1[0] / 2;
            tempstate[1] = _state.lp + deriv1[1] / 2;

            np = input - 2 * R * _state.bp;
            hp = np - _state.lp;
            deriv2[0] = h * TWOPI * w_c * trig.tanh_LUT(hp); // is this the right move?
            deriv2[1] = h * TWOPI * w_c * trig.tanh_LUT(tempstate[0]);

            tempstate[0] = _state.bp + deriv2[0] / 2;
            tempstate[1] = _state.lp + deriv2[1] / 2;

            np = input - 2 * R * _state.bp;
            hp = np - _state.lp;
            deriv3[0] = h * TWOPI * w_c * trig.tanh_LUT(hp);
            deriv3[1] = h * TWOPI * w_c * trig.tanh_LUT(tempstate[0]);

            tempstate[0] = _state.bp + deriv3[0];
            tempstate[1] = _state.lp + deriv3[1];

            np = input - 2 * R * _state.bp;
            hp = np - _state.lp;
            deriv4[0] = h * TWOPI * w_c * trig.tanh_LUT(hp);
            deriv4[1] = h * TWOPI * w_c * trig.tanh_LUT(tempstate[0]);

            _state.bp += (1.f/6.f) * (deriv1[0] + 2 * deriv2[0] + 2 * deriv3[0] + deriv4[0]);
            _state.lp += (1.f/6.f) * (deriv1[1] + 2 * deriv2[1] + 2 * deriv3[1] + deriv4[1]);
        }
        
        _outputs.bp = _state.bp;
        _outputs.lp = _state.lp;
        _outputs.hp = hp;
        _outputs.np = np;
    }
private:
    int _oversample_factor;
    float h, w_c, R, resonance;
    // k_1 through k_4. for each, [0] is bp, [1] is lp.
    float deriv1[2], deriv2[2], deriv3[2], deriv4[2]; 
};

//==============================================================================
class slewlim
{
public:
    slewlim()
        : _vOut(0.f)
    {
        setSampleRate(44100.f);
    }
    slewlim(float sample_rate)
        : _vOut(0.f)
    {
        setSampleRate(sample_rate);
    }
    ~slewlim()
    {
    }
    //============================================================
    void setSampleRate(float sample_rate)
    {
        this->sampleRate = sample_rate;
        this->T = 1.f / sample_rate;
    }
    //============================================================
    // immediate change
    void setRise(float rise)
    {
        this->rise = rise;
        this->riseInc = (this->T * 1000.f) / (this->rise);
    }
    // change over block size
    void setRise()
    {
        if (_riseTarget != rise)
        {
            this->rise += (_riseTarget - rise) * _oneOverBlockSize;
            this->riseInc = (this->T * 1000.f) / (rise);
        }
    }
    void setRise(float riseTarget, float oneOverBlockSize)
    {
        this->rise += (riseTarget - this->rise) * oneOverBlockSize;
        this->riseInc = (this->T * 1000.f) / (this->rise);
    }
    // immediate change
    void setFall(float fall)
    {
        this->fall = fall;
        this->fallInc = (this->T * 1000.f) / (this->fall);
    }
    // change over block size
    void setFall()
    {
        if (_fallTarget != fall)
        {
            this->fall += (_fallTarget - fall) * _oneOverBlockSize;
            this->fallInc = (this->T * 1000.f) / (fall);
        }
    }
    void setFall(float fallTarget, float oneOverBlockSize)
    {
        this->fall += (fallTarget - this->fall) * oneOverBlockSize;
        this->fallInc = (this->T * 1000.f) / (this->fall);
    }
    //============================================================
    float ASR(float gate)
    {
        using namespace nkvdu_memoryless;
        
        setRise();
        setFall();
        
        if (_vOut < gate)
        {
            _vOut += riseInc;
            _vOut = clamp_high<float>(_vOut, gate);
        }
        else if (_vOut > gate)
        {
            _vOut -= fallInc;
            _vOut = clamp_low<float>(_vOut, gate);
        }
        
        return _vOut;
    }

    float _riseTarget, _fallTarget;
    float _oneOverBlockSize;
private:
    float sampleRate, T;
    
    // 'Inc' variables tell change per sample.
    float rise, riseInc, fall, fallInc, _vOut;
};

class dcBlock   :   public filter_abstract
{
public:
    dcBlock()
    {
        setSampleRate(44100);
        clear();
        R = 0.995;  //fixed R
    }
    dcBlock(float sample_rate)
    {
        setSampleRate(sample_rate);
        clear();
        R = 0.995;  //fixed R, should base on sample rate instead
    }
    void clear()
    {
        xz1 = yz1 = 0.f;
    }
    void updateCutoff(float cutoff_target, float oneOverBlockSize)
    {
    }
    void updateResonance(float res_target, float oneOverBlockSize)
    {
    }
    void updateR(float R_target, float oneOverBlockSize)
    {
        R += (R_target - R) * oneOverBlockSize;
    }
    float filter(float x)
    {
        yz1 = x - xz1 + R * yz1;
        xz1 = x;
        return yz1;
    }
private:
    float R, xz1, yz1;
};

//============================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
class CTPTMoogFilterStage
{
public:
    CTPTMoogFilterStage()
    {
    }
    ~CTPTMoogFilterStage()
    {
    }
protected:
    // controls
    float G;
    float scalar;
    float sampleRate;
    float z1;

    // cutoff
    // scalar value
    // fs
    // z-1 storage location
public:
    inline void initialize(float newSampleRate)
    {
        // save
        sampleRate = newSampleRate;
        z1 = 0;
    }
    void setFc(float fc)
    {
        // prewarp the cutoff- these are bilinear-transform filters
        float wd = 2 * PI * fc;
        float T  = 1 / sampleRate;
        float wa = (2 / T) * tan(wd * T / 2);
        float g  = wa * T / 2;
        // calculate big G value; see Zavalishin p46 The Art of VA Design
        G = g / (1.0 + g);
    }
    
    float doFilterStage(float xn)
    {
        float v = (xn - z1) * G;
        float out = v + z1;
        z1 = out + v;
        return out;
    }
    float getSampleRate()
    {
        return sampleRate;
    }
    float getStorageRegisterValue()
    {
        return z1;
    }
};
class CTPTMoogLadderFilter
{
public:
    CTPTMoogLadderFilter()
    {
    }
    ~CTPTMoogLadderFilter()
    {
    }
protected:
    CTPTMoogFilterStage filter1;
    CTPTMoogFilterStage filter2;
    CTPTMoogFilterStage filter3;
    CTPTMoogFilterStage filter4;
    float k; // Q control
    float fc; // fc control
public:
    inline void initialize(float newSampleRate)
    {
        filter1.initialize(newSampleRate);
        filter2.initialize(newSampleRate);
        filter3.initialize(newSampleRate);
        filter4.initialize(newSampleRate);
    }
    inline void calculateTPTCoeffs(float cutoff, float Q)
    {
        // 4 sync-tuned filters
        filter1.setFc(cutoff);
        filter2.setFc(cutoff);
        filter3.setFc(cutoff);
        filter4.setFc(cutoff);
        
        // NOTE: Q is limited to 20 on the UI to prevent blowing up //Q=0.707->25==>k=0->4
        k = 4.0*(Q - 0.707)/(25.0 - 0.707);
        // ours
        fc = cutoff;
    }
    float doTPTMoogLPF(float xn)
    {
        // calculate g
        float wd = 2 * PI * fc;
        float T  = 1 / (float)filter1.getSampleRate();
        float wa = (2 / T) * tan(wd * T / 2);
        float g  = wa * T / 2;
        float G = g * g * g * g;
        float S = g * g * g * filter1.getStorageRegisterValue() +
        g * g * filter2.getStorageRegisterValue() +
        g * filter3.getStorageRegisterValue() +
        filter4.getStorageRegisterValue();
        //uis input to filters, straight from book
        float u = (xn - k * S) / (1 + k * G);
        // four cascades using nested functions
        float filterOut = filter4.doFilterStage(filter3.doFilterStage (filter2.doFilterStage(filter1.doFilterStage(u))));
        // output
        return filterOut;
    }
};

class tvap  :   public filter_abstract
{
public:
    tvap() {
    }
    ~tvap() {
    }
    //==============================================================================
    virtual void setSampleRate(float sample_rate)
    {
        this->sampleRate = sample_rate;
        this->T = 1.f / sample_rate;
        lp.setSampleRate(sample_rate);
        lp.updateCutoff(sample_rate * 0.125f, 1.f);
    }
    
    void clear(){
        sp->z1 = sp->z2 = sp->fb_proc = 0.f;
    };
    void updateCutoff(float cutoff_target, float oneOverBlockSize) {
        if (cutoff_target < 0) cutoff_target = 0;
        f_pi += (cutoff_target - f_pi) * oneOverBlockSize;
        calc_b1();
    }
    void updateResonance(float res_target, float oneOverBlockSize) {
        //if (f_b <= 0) f_b = 0.0000000001;
        f_b += (res_target - f_b) * oneOverBlockSize;
        f_b_to_b0();
    }
// function aliases just to have more meaningful names
    void update_f_pi (float f_pi_target, float oneOverBlockSize)
    {   updateCutoff(f_pi_target, oneOverBlockSize); }
    void update_f_b(float f_b_target, float oneOverBlockSize)
    {   updateResonance(f_b_target, oneOverBlockSize); }
    
    void calc_b1(void) {
        float d = -1 * cosf((2.f * PI * f_pi) / sampleRate);
        float c = (tanf(PI * f_b / sampleRate) - 1.f) / (tanf(PI * f_b / sampleRate) + 1.f);
        float r1 = acosf(-1.f * c);
        float r2 = acosf(-1.f * d);
        b1 = cosf(r1) * (1.f + cosf(r2));
    }
    void f_b_to_b0(void) {
        float c = (tanf(PI * f_b / sampleRate) - 1.f) / (tanf(PI * f_b / sampleRate) + 1.f);
        float r1 = acosf(-1.f * c);
        b0 = cosf(r1);
    }

    float f_pi2r2(float _f_pi)
    {
        float d = -1 * cosf((2.f * PI * _f_pi) / sampleRate);
        float r2 = acosf(-d);
        return r2;
    }
    float f_b2r1(float _f_b)
    {
        float tmp = tanf(PI * _f_b / sampleRate);
        float c = (tmp - 1.f) / (tmp + 1.f);
        float r1 = acosf(-c);
        return r1;
    }

    float filter(float x_n) {
        /* float _y1 = state.y1;
        float _y2 = state.y2;
        float _x1 = state.x1;
        float _x2 = state.x2;
        float y_n = b0 * x_n - b1 * _x1 + _x2 + b1 * _y1 - b0 * _y2;
        state.y2 = _y1;
        state.y1 = y_n;
        state.x2 = _x1;
        state.x1 = x_n;
        return y_n; */
        float _r1, _r2, _cr1, _cr2, _sr1, _sr2;
        float tmp [3];
        _r1 = f_b2r1(f_b);
        _r2 = f_pi2r2(f_pi);

        _cr1 = cos(_r1);
        _sr1 = sin(_r1);
        _cr2 = cos(_r2);
        _sr2 = sin(_r2);
        //tmp[0] = _x_n;
        tmp[1] = _cr2 * z1 - _sr2 * state.z2;
        //tmp[2] = _sr2 * _z1 + _cr2 * _z2;

        tmp[0] = _cr1 * x_n - _sr1 * tmp[1];
        tmp[1] = _sr1 * x_n + _cr1 * tmp[1];
        tmp[2] = _sr2 * z1 + _cr2 * state.z2;

        state.z1 = tmp[1];
        state.z2 = tmp[2];

        return tmp[0];
    }
    
    // should be in memoryless but got linker error
    inline float unboundSat2(float x)
    {
        float num = 2.f * x;
        float denom = 1.f + sqrt(1.f + fabs(4.f * x));
        return num / denom;
    }
    
    float filter_fbmod(float x_n, float fb_f_pi, float fb_f_b)
    {
        float _f_pi_n, _f_b_n;
        _f_pi_n = f_pi;
        _f_b_n = f_b;

        fb_f_pi *= f_pi;
        fb_f_pi *= 1000.f;
        fb_f_b *= f_pi;
        fb_f_b *= 1000.f;
        
        _f_pi_n += unboundSat2(state.fb_proc * fb_f_pi)* 10.f;
        _f_b_n  += unboundSat2(state.fb_proc * fb_f_b) * 10.f;

        // if (_f_b_n < 0.f) _f_b_n = 0.f;
        if (_f_b_n < 0.f) 
            _f_b_n = expf(_f_b_n);
        else    
            _f_b_n += 1.f;

        float highLimit = (sampleRate * 0.495f);
        if (_f_pi_n >= highLimit) _f_pi_n = highLimit;
        if (_f_b_n >= highLimit) _f_b_n = highLimit;

        float _r1, _r2, _cr1, _cr2, _sr1, _sr2;
        float tmp[3];
        _r1 = f_b2r1(_f_b_n);
        _r2 = f_pi2r2(_f_pi_n);

        _cr1 = cos(_r1);
        _sr1 = sin(_r1);
        _cr2 = cos(_r2);
        _sr2 = sin(_r2);
        //tmp[0] = _x_n;
        tmp[1] = _cr2 * state.z1 - _sr2 * state.z2;
        //tmp[2] = _sr2 * _z1 + _cr2 * _z2;

        tmp[0] = _cr1 * x_n - _sr1 * tmp[1];
        tmp[1] = _sr1 * x_n + _cr1 * tmp[1];
        tmp[2] = _sr2 * state.z1 + _cr2 * state.z2;

        state.z1 = tmp[1];
        state.z2 = tmp[2];

        
        float fb_filt = dcFilt.filter(tmp[0]);   // used to feed output to modulate control inputs
        fb_filt = lp.tpt_lp(fb_filt);
        //fb_filt = 0.5f * (fb_filt + _fb_filt_z1);
        state.fb_proc = fb_filt;

        return tmp[0];
    }
    
protected:
    typedef struct tvapstate {
        //float x1, x2, y1, y2;
        float z1, z2;
        float fb_proc;  // processed fed back output sample
    } _tvapstate;
    _tvapstate state = {.z1 = 0.f, .z2 = 0.f, .fb_proc =0.f};
    _tvapstate *sp = &state;

    dcBlock dcFilt;
    onePole lp;
private:
    float f_pi, f_b;
    float b0, b1;
};

}
