/* Copyright (c) 2017 Kevin Nelson
 * Based on line~, which is Copyright (c) 1997-1999 Miller Puckette.
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

#include "m_pd.h"  
#include "math.h"

/* -------------------------- line~ ------------------------------ */
static t_class *expline_tilde_class;

typedef struct _expline
{
    t_object x_obj;
    t_sample x_target; /* output target value of ramp */
    t_sample x_value; /* current output value of ramp at block-borders */

    double x_scale_value; /* current value of 0-1 ramp at block-borders */
    double x_overshoot_target;
    double x_overshoot_ratio;
    double x_attack_coef;

    t_float c_last_samplerate;
    t_float c_last_dspticksize;
    t_float x_dspticktomsec;
    t_float x_samplespermsec;
    int x_needs_overshoot_recalc;

    t_sample x_valoffset;
    t_sample x_valmult;

    t_float x_inlet_ramptime;
    t_float x_inlet_ramptime_was;
    t_float x_inlet_overshoot;
    t_float x_inlet_overshoot_was;
    int x_ticksleft;
    int x_retarget;
} t_expline;

/* M_E is `e`.. and a double */

static void recalc_attack(t_expline *x)
{
    post("re-calculating attack");
    post("sample rate old: %f new: %f", x->c_last_samplerate, x->x_samplespermsec);
    post("dsptick old: %f new: %f", x->c_last_dspticksize, x->x_dspticktomsec);
    double attack_samples = x->x_samplespermsec * x->x_inlet_ramptime_was;
    x->x_attack_coef = 1.0 - exp(log(x->x_overshoot_ratio) / attack_samples);
    post("attack samples: %f", attack_samples);
    post("attack coef: %f", x->x_attack_coef);
    x->x_needs_overshoot_recalc = 0;
}

static t_int *expline_tilde_perform(t_int *w)
{
    t_expline *x = (t_expline *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    t_sample f = x->x_value;
    double sv = x->x_scale_value;

    if (PD_BIGORSMALL(f))
            x->x_value = f = 0;
    if (x->x_retarget)
    {
        int nticks = x->x_inlet_ramptime_was * x->x_dspticktomsec;
        if (!nticks) nticks = 1;
        x->x_ticksleft = nticks;
        recalc_attack(x);
        x->x_retarget = 0;
    }
    if (x->x_ticksleft)
    {
        if (x->x_needs_overshoot_recalc == 1 ||
                (x->c_last_samplerate != x->x_samplespermsec) ||
                (x->c_last_dspticksize != x->x_dspticktomsec)) {
            recalc_attack(x);
        }
        t_sample g = x->x_value;
        while (n--) {
            *out++ = g;
            sv += x->x_attack_coef * (x->x_overshoot_target - sv);
            if (sv > 1.0) sv = 1.0;
            g = x->x_valoffset + (x->x_valmult * sv);
        }
        x->x_value = g;
        x->x_scale_value = sv;
        x->x_ticksleft--;
    }
    else
    {
        t_sample g = x->x_value = x->x_target;
        while (n--)
            *out++ = g;
    }

    /* Cache the dsp block size / sample rate to determine if we need
     * to recalculate the attack coeffient
     */
    x->c_last_samplerate = x->x_samplespermsec;
    x->c_last_dspticksize = x->x_dspticktomsec;
    return (w+4);
}

/* TB: vectorized version */
static t_int *expline_tilde_perf8(t_int *w)
{
    t_expline *x = (t_expline *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    t_sample f = x->x_value;
    double sv = x->x_scale_value;

    if (PD_BIGORSMALL(f))
        x->x_value = f = 0;
    if (x->x_retarget)
    {
        int nticks = x->x_inlet_ramptime_was * x->x_dspticktomsec;
        if (!nticks) nticks = 1;
        x->x_ticksleft = nticks;
        recalc_attack(x);
        x->x_retarget = 0;
    }
    if (x->x_ticksleft)
    {
        if (x->x_needs_overshoot_recalc == 1 ||
                (x->c_last_samplerate != x->x_samplespermsec) ||
                (x->c_last_dspticksize != x->x_dspticktomsec)) {
            recalc_attack(x);
        }
        t_sample g = x->x_value;
        while (n--) {
            *out++ = g;
            sv += x->x_attack_coef * (x->x_overshoot_target - sv);
            if (sv > 1.0) sv = 1.0;
            g = x->x_valoffset + (x->x_valmult * sv);
        }
        x->x_value = g;
        x->x_scale_value = sv;
        x->x_ticksleft--;
    }
    else
    {
        t_sample g = x->x_value = x->x_target;
        for (; n; n -= 8, out += 8)
        {
            out[0] = g; out[1] = g; out[2] = g; out[3] = g; 
            out[4] = g; out[5] = g; out[6] = g; out[7] = g;
        }
    }

    /* Cache the dsp block size / sample rate to determine if we need
     * to recalculate the attack coeffient
     */
    x->c_last_samplerate = x->x_samplespermsec;
    x->c_last_dspticksize = x->x_dspticktomsec;
    return (w+4);
}

static void set_overshoot(t_expline *x, t_float overshoot)
{
    x->x_inlet_overshoot_was = overshoot;
    x->x_overshoot_target = 1.0 + overshoot;
    x->x_overshoot_ratio = overshoot / x->x_overshoot_target;
    post("overshoot target: %f", x->x_overshoot_target);
    x->x_needs_overshoot_recalc = 1;
}

static void expline_tilde_float(t_expline *x, t_float f)
{
    post("previous attack coef: %f", x->x_attack_coef);
    if (x->x_inlet_overshoot != x->x_inlet_overshoot_was)
    {
        set_overshoot(x, x->x_inlet_overshoot);
    }

    if (x->x_inlet_ramptime <= 0)
    {
        /* handle the "set to this value now" case */
        x->x_target = x->x_value = f;
        x->x_ticksleft = x->x_retarget = x->x_scale_value = 0;
    }
    else
    {
        /* ramp to a new target */
        x->x_target = f;
        x->x_valoffset = x->x_value;
        x->x_valmult = f - x->x_value;
        x->x_retarget = 1;
        x->x_scale_value = 0.0;
        x->x_inlet_ramptime_was = x->x_inlet_ramptime;
        x->x_inlet_ramptime = 0;
        post("ramp time: %f", x->x_inlet_ramptime_was);
    }
}

static void expline_tilde_stop(t_expline *x)
{
    x->x_target = x->x_value;
    x->x_ticksleft = x->x_retarget = x->x_scale_value = 0;
}

static void expline_tilde_dsp(t_expline *x, t_signal **sp)
{
    if(sp[0]->s_n&7)
        dsp_add(expline_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
    else
        dsp_add(expline_tilde_perf8, 3, x, sp[0]->s_vec, sp[0]->s_n);

    x->x_samplespermsec = sp[0]->s_sr / 1000;
    x->x_dspticktomsec = x->x_samplespermsec / sp[0]->s_n;
}

static void *expline_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_expline *x = (t_expline *)pd_new(expline_tilde_class);
    outlet_new(&x->x_obj, gensym("signal"));
    floatinlet_new(&x->x_obj, &x->x_inlet_ramptime);
    floatinlet_new(&x->x_obj, &x->x_inlet_overshoot);
    x->x_ticksleft = x->x_retarget = 0;
    x->x_scale_value = x->x_value = x->x_target = x->x_inlet_ramptime = x->x_inlet_ramptime_was = x->x_inlet_overshoot = x->x_inlet_overshoot_was = 0;

    t_float overshoot = 0.1;

    switch(argc) {
        case 1:
            overshoot = atom_getfloat(argv);
            break;
        default:
            break;
    }

    set_overshoot(x, overshoot);
    return (x);
}

void expline_tilde_setup(void)
{
    expline_tilde_class = class_new(gensym("expline~"), (t_newmethod)expline_tilde_new, 0,
        sizeof(t_expline), CLASS_DEFAULT, A_GIMME, 0);
    class_addfloat(expline_tilde_class, (t_method)expline_tilde_float);
    class_addmethod(expline_tilde_class, (t_method)expline_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(expline_tilde_class, (t_method)expline_tilde_stop,
        gensym("stop"), 0);
}
