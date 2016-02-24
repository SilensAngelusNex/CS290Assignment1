===============
CS290Assignment1
================
Aaron Liberatore
Weston Carvalhosd



## Near vs Far ##



## Ellipsoid ##

#### Calculations #####

a = 10 b = 7 c =7
x' = ax
y' = by
z' = cz
(z')2 / 100 + (y')2 / 49 + (z')2 / 49 = 1

100x2 / 100 + 49y2 / 49 = 1

a^2 - c^2 = b^2
100 - c^2 = 49
100-49 = c^2
51 = 7.14

Foci: (7.14,0,0) and (-7.14,0,0)

#### Results #####
When the paths are drawn with the source and the receiver set to the two foci of the ellipse, it appears as though there are 1 order bounces from everyface to the other foci. It becomes incredibly hard to actually see anything because of all the paths. Below is a clearer view from the source itself.

![Source and receiver at foci](./points_at_foci.png "Paths are pink")

Despite there being a path from source to receiver off everyface of the ellipsoid when each is at the foci, if you move the source/receiver over by merely a couple meters there are only a few paths.

![Source and receiver at foci](./off_foci.png "Paths are pink")
