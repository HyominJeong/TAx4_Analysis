
# coding: utf-8

# In[2]:


import numpy as np


# In[4]:


g = 9.8 # m/s**2


# In[13]:


def throw(theta_rad, v0 = 10):
    v0_x = v0 * np.cos(theta_rad)
    v0_y = v0 * np.sin(theta_rad)
    t = 2. * v0_y / g
    x = v0_x * t
    return x


# In[29]:


x = 0. 


# In[31]:


while (x < 4.95) or (x > 5.05):
    # When you want to input theta
    theta = float(input("Input theta (degree unit, 0 to 90)"))
    #break
    theta_rad = theta * np.pi / 180
    x = throw(theta_rad)
    if x > 5.05:
        print("Too long")
    elif x < 4.95:
        print("Too short")
    else:
        print("Goal, theta is %2.2f (deg) and distance is %2.2f (m)" % (theta, x))

