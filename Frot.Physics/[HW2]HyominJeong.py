
# coding: utf-8

# In[1]:


import numpy as np


# In[9]:


# Define input data set
# distance(Mpc unit), velocity(km/s) set
dis_vel_set = np.array([[15, 1100],

[97, 6700],

[32, 2400],

[145, 10700],

[50, 3100],

[122, 9900],

[58, 4300],

[91, 5300],

[120, 9000],

[93, 7500],

[158, 8900],

[64, 5300],

[145, 9600],

[61, 3300],

[103, 5100],

[46, 3600],

[34, 1800],

[185, 9500],

[20, 1200]])


# In[15]:


print(dis_vel_set.transpose()[0])


# In[16]:


distance_set = dis_vel_set.transpose()[0]
velocity_set = dis_vel_set.transpose()[1]
print(distance_set)
print(velocity_set)


# In[114]:


# by Hubble law, vel = H * dis
H_init = np.random.random() * 100 # trial
print("Initial Hubble constant is %.2f" % H_init)

H_try = H_init


# In[69]:


def vel_error_sq(distance, velocity):
    return (velocity - H_try * distance)**2


# In[70]:


def vel_error(distance, velocity):
    return H_try * distance - velocity


# In[121]:


for i in range(100000):
    alpha = 0.0001 / (i+1) # Learning rate, decrease with the number of learning
    # Learn
    for dis_vel in dis_vel_set:
        #print(dis_vel, H_try)
        error = vel_error(dis_vel[0], dis_vel[1])
        #print(error)
        H_try = H_try - alpha * error * dis_vel[0]
    phi = 1.* sum(vel_error_sq(distance_set, velocity_set)) / len(distance_set)
    if (i % 10000) == 0:
        print("Hubble constant is %5.2f, phi is %2.2f, %6i/%i Learned" % (H_try, phi, i, 100000))


# In[117]:


print("Hubble constant is %.2f" % H_try)

