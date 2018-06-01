# A shallow water equation model

------

This is a shallow water equation model written in python. I wrote this model only for practice of numerical simulation and python coding following a youtube video by Arun Prasaad Gunasekaran (https://www.youtube.com/watch?v=7_kcM8Yfv3A). So the framework and even code are nearly identical to his model which is written using Fortran. The only difference is the use of 1.0/(2.0*dx) instead of 1.0/2.0*dx when calcualte horizontal gradient (this might be a typo remaining in Arun's code).

The basic feature of this model including:

> * Non-conservative form shallow water equations (referred to https://en.wikipedia.org/wiki/Shallow_water_equations)
> * Coriolis, frictional and viscous force are included
> * Arakawa-A grid using center finite difference scheme
> * Forth-order Runge-Kutta methods for time integration

![shallowWaterEquation](https://wikimedia.org/api/rest_v1/media/math/render/svg/87347a5f2337ec3c8042bac5274e5837d1047f83)

Because I'm a newbee in numerical computing and python coding, massive problems remained in the code. Any suggestion on the model, e.g., how to improve the calculation efficient, or how to make the model more realistic, are welcome.

Special thanks to Arun Prasaad Gun for his youtube video on coding. It's really helpful for me.
