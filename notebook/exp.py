
# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab inline')
pylab.rcParams['figure.figsize'] = (18, 6)

import numpy as np


# In[2]:

sample_rate = 48000
tau = 2.0 * np.pi


# In[3]:

def exp_ramp(current, target, attack, linearity):
    flip = False
    if current > target:
        flip = True
        current, target = target, current
    out = current
    out_arr = [out]
    sample_count = 0
    for i in range(0, 500):
        out += attack * ((target/linearity) - out)
        if flip:
            out_arr.append(target - out + current)
        else:
            out_arr.append(out)
        if out > target:
            out = target
        else:
            sample_count += 1
    print("reach %s in %s samples" % (target, sample_count))
    return np.array(out_arr)


# In[4]:

plot(exp_ramp(0.0, 1.0, 0.006, 0.63))


# In[5]:

plot(exp_ramp(0.0, 1.0, 0.003, 0.63))


# In[6]:

plot(exp_ramp(0.0, 1.0, 0.003, 0.7))


# In[7]:

plot(exp_ramp(0.0, 1.0, 0.004, 0.54))


# In[8]:

plot(exp_ramp(0.0, 1.0, 0.008, 0.85))


# In[9]:

plot(exp_ramp(0.0, 1.0, 0.008, 0.85))


# In[10]:

plot(exp_ramp(0.5, 1.0, 0.008, 0.85))


# In[11]:

plot(exp_ramp(0.0, 0.5, 0.005, 0.5))


# In[12]:

plot(exp_ramp(0.0, 1.0, 0.01, 0.875))


# In[ ]:




# In[13]:

def exp_ramp2(start, end, attack, linearity):
    # always start at 0, so the curve of the ramp stays consistent
    offset = start
    out_mult = end - start
    
    out = 0.0
    out_arr = [out]
    sample_count = 0
    for i in range(0, 500):
        out += attack * ((1.0/linearity) - out)
        out_arr.append(offset + (out_mult * out))
        if out > 1.0:
            out = 1.0
        else:
            sample_count += 1
    print("reach %s in %s samples" % (1.0, sample_count))
    return np.array(out_arr)


# In[14]:

plot(exp_ramp2(0.0, 1.0, 0.005, 0.5))


# In[15]:

plot(exp_ramp2(0.5, 1.0, 0.005, 0.5))


# In[16]:

plot(exp_ramp2(1.0, 0.1, 0.005, 0.5))


# In[17]:

plot(exp_ramp2(1.0, 0.1, 0.0088, 0.88))


# In[18]:

sample_rate = 44100
attack_ms = 100
attack_samples = sample_rate / 1000.0 * attack_ms
print("attack samples: %s" % attack_samples)


# In[19]:

def attack_coef_from_ms(linearity_coef, attack_ms, sample_rate):
    attack_samples = sample_rate / 1000.0 * attack_ms
    coef = np.e * linearity_coef * linearity_coef / attack_samples
    return coef
print(attack_coef_from_ms(0.88, 100, 44100))


# In[20]:

def exp_ramp_ms(start, end, attack_ms, linearity, sample_rate):
    attack_samples = floor(sample_rate / 1000.0 * attack_ms)
    attack_coef = np.e * linearity * linearity / attack_samples
    print("attack coef: %s" % attack_coef)
    # always start at 0, so the curve of the ramp stays consistent
    offset = start
    out_mult = end - start
    
    out = 0.0
    out_arr = [out]
    sample_count = 0
    while 1:
        out += attack_coef * ((1.0/linearity) - out)
        out_arr.append(offset + (out_mult * out))
        if out > 1.0:
            out = 1.0
            break
        else:
            sample_count += 1
    print("reached target in %s samples" % sample_count)
    return np.array(out_arr)


# In[21]:

plot(exp_ramp_ms(0.0, 1.0, 100, 0.88, 44100))


# In[22]:

plot(exp_ramp_ms(0.0, 1.0, 30, 0.88, 44100))


# In[23]:

plot(exp_ramp_ms(0.0, 1.0, 30, 0.8, 44100))


# In[24]:

plot(exp_ramp_ms(0.0, 1.0, 30, 0.7, 44100))


# In[25]:

def exp_ramp_ms_hard(start, end, attack_ms, linearity, sample_rate):
    attack_samples = int(floor(sample_rate / 1000.0 * attack_ms))
    attack_coef = np.e * linearity * linearity / attack_samples
    print("attack coef: %s" % attack_coef)
    # always start at 0, so the curve of the ramp stays consistent
    offset = start
    out_mult = end - start
    
    out = 0.0
    out_arr = [out]
    sample_count = 0
    for i in range(0, attack_samples):
        out += attack_coef * ((1.0/linearity) - out)
        out_arr.append(offset + (out_mult * out))
        if out > 1.0:
            out = 1.0
            break
        else:
            sample_count += 1
    print("reached %s in %s samples" % (out, sample_count))
    return np.array(out_arr)


# In[26]:

plot(exp_ramp_ms_hard(0.0, 1.0, 30, 0.88, 44100))


# In[27]:

plot(exp_ramp_ms_hard(0.0, 1.0, 500, 0.44, 44100))


# In[28]:

plot(exp_ramp_ms_hard(0.0, 1.0, 500, 0.63, 44100))


# In[31]:

def plot_lin(start, end, time_ms, sample_rate):
    samples = sample_rate / 1000 * time_ms
    increment = float(end - start) / samples
    out = start
    out_arr = [out]
    for i in range(0, samples):
        out += increment
        out_arr.append(out)
    return np.array(out_arr)


# In[32]:

plot(exp_ramp_ms_hard(0.0, 1.0, 30, 0.88, 44100))
plot(plot_lin(0.0, 1.0, 30, 44100))


# In[33]:

plot(exp_ramp_ms(0.0, 1.0, 30, 0.8, 44100))
plot(plot_lin(0.0, 1.0, 30, 44100))


# In[34]:

plot(exp_ramp_ms(0.0, 1.0, 30, 0.63, 44100))
plot(plot_lin(0.0, 1.0, 30, 44100))


# In[35]:

plot(exp_ramp_ms(0.0, 1.0, 100, 0.63, 44100))
plot(plot_lin(0.0, 1.0, 100, 44100))


# In[36]:

plot(exp_ramp_ms(0.0, 1.0, 500, 0.63, 44100))
plot(plot_lin(0.0, 1.0, 500, 44100))


# In[37]:

plot(exp_ramp_ms(0.0, 1.0, 500, 0.75, 44100))
plot(plot_lin(0.0, 1.0, 500, 44100))


# In[38]:

plot(exp_ramp_ms(0.0, 1.0, 500, 0.44, 44100))
plot(plot_lin(0.0, 1.0, 500, 44100))


# In[41]:

plot(exp_ramp_ms(0.0, 1.0, 500, 0.525, 44100))
plot(plot_lin(0.0, 1.0, 500, 44100))


# In[61]:

plot(exp_ramp_ms(0.0, 1.0, 500, 0.7, 48000))
plot(plot_lin(0.0, 1.0, 500, 48000))


# In[65]:

plot(exp_ramp_ms(0.0, 1.0, 10, 0.55, 48000))
plot(plot_lin(0.0, 1.0, 10, 48000))


# In[66]:

plot(exp_ramp_ms(67, 73, 10, 0.55, 48000))
plot(plot_lin(67, 73, 10, 48000))


# In[ ]:



