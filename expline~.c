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
    t_sample x_target; /* target value of ramp */
    t_sample x_value; /* current value of ramp at block-borders */
    double x_overshoot_mult;
    double x_attack_coef;
    t_sample x_biginc;
    t_sample x_inc;
    t_float x_1overn;
    t_float x_dspticktomsec;
    t_float x_inlet_ramptime;
    t_float x_inlet_ramptime_was;
    t_float x_inlet_overshoot;
    t_float x_inlet_overshoot_was;
    int x_ticksleft;
    int x_retarget;
} t_expline;

/* M_E is `e`.. and a double */

static t_int *expline_tilde_perform(t_int *w)
{
    t_expline *x = (t_expline *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    t_sample f = x->x_value;

    if (PD_BIGORSMALL(f))
            x->x_value = f = 0;
    if (x->x_retarget)
    {
        int nticks = x->x_inlet_ramptime_was * x->x_dspticktomsec;
        if (!nticks) nticks = 1;
        x->x_ticksleft = nticks;
        x->x_biginc = (x->x_target - x->x_value)/(t_float)nticks;
        x->x_inc = x->x_1overn * x->x_biginc;
        x->x_retarget = 0;
    }
    if (x->x_ticksleft)
    {
        t_sample f = x->x_value;
        while (n--) *out++ = f, f += x->x_inc;
        x->x_value += x->x_biginc;
        x->x_ticksleft--;
    }
    else
    {
        t_sample g = x->x_value = x->x_target;
        while (n--)
            *out++ = g;
    }
    return (w+4);
}

/* TB: vectorized version */
static t_int *expline_tilde_perf8(t_int *w)
{
    t_expline *x = (t_expline *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    t_sample f = x->x_value;

    if (PD_BIGORSMALL(f))
        x->x_value = f = 0;
    if (x->x_retarget)
    {
        int nticks = x->x_inlet_ramptime_was * x->x_dspticktomsec;
        if (!nticks) nticks = 1;
        x->x_ticksleft = nticks;
        x->x_biginc = (x->x_target - x->x_value)/(t_sample)nticks;
        x->x_inc = x->x_1overn * x->x_biginc;
        x->x_retarget = 0;
    }
    if (x->x_ticksleft)
    {
        t_sample f = x->x_value;
        while (n--) *out++ = f, f += x->x_inc;
        x->x_value += x->x_biginc;
        x->x_ticksleft--;
    }
    else
    {
        t_sample f = x->x_value = x->x_target;
        for (; n; n -= 8, out += 8)
        {
            out[0] = f; out[1] = f; out[2] = f; out[3] = f; 
            out[4] = f; out[5] = f; out[6] = f; out[7] = f;
        }
    }
    return (w+4);
}

static void expline_tilde_float(t_expline *x, t_float f)
{
    if (x->x_inlet_overshoot != x->x_inlet_overshoot_was)
    {
        x->x_inlet_overshoot_was = x->x_inlet_overshoot;
        double overshoot_dest = 1.0 + x->x_inlet_overshoot;
        x->x_overshoot_mult = M_E * overshoot_dest * overshoot_dest;
    }

    if (x->x_inlet_ramptime <= 0)
    {
        x->x_target = x->x_value = f;
        x->x_ticksleft = x->x_retarget = 0;
    }
    else
    {
        x->x_target = f;
        x->x_retarget = 1;
        x->x_inlet_ramptime_was = x->x_inlet_ramptime;
        x->x_inlet_ramptime = 0;
    }
}

static void expline_tilde_stop(t_expline *x)
{
    x->x_target = x->x_value;
    x->x_ticksleft = x->x_retarget = 0;
}

static void expline_tilde_dsp(t_expline *x, t_signal **sp)
{
    if(sp[0]->s_n&7)
        dsp_add(expline_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
    else
        dsp_add(expline_tilde_perf8, 3, x, sp[0]->s_vec, sp[0]->s_n);
    x->x_1overn = 1./sp[0]->s_n;
    x->x_dspticktomsec = sp[0]->s_sr / (1000 * sp[0]->s_n);
}

static void *expline_tilde_new(t_floatarg overshoot)
{
    t_expline *x = (t_expline *)pd_new(expline_tilde_class);
    x->x_overshoot = overshoot;
    outlet_new(&x->x_obj, gensym("signal"));
    floatinlet_new(&x->x_obj, &x->x_inlet_ramptime);
    floatinlet_new(&x->x_obj, &x->x_inlet_overshoot);
    x->x_ticksleft = x->x_retarget = 0;
    x->x_value = x->x_target = x->x_inlet_ramptime = x->x_inlet_ramptime_was = x->x_inlet_overshoot = x->x_inlet_overshoot_was = 0;
    return (x);
}

void expline_tilde_setup(void)
{
    expline_tilde_class = class_new(gensym("expline~"), expline_tilde_new, 0,
        sizeof(t_expline), CLASS_DEFAULT, A_DEFFLOAT, 0.1);
    class_addfloat(expline_tilde_class, (t_method)expline_tilde_float);
    class_addmethod(expline_tilde_class, (t_method)expline_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(expline_tilde_class, (t_method)expline_tilde_stop,
        gensym("stop"), 0);
}
