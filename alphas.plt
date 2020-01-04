set fit errorvariables

m0=0.0; c0=4.62497

f1(x)=-m1*x+c1
fit f1(x) '$0.ISF' using 1:3 via m1, c1

f2(x)=-m2*x+c2
fit f2(x) '$0.ISF' using 1:4 via m2, c2

f3(x)=-m3*x+c3
fit f3(x) '$0.ISF' using 1:5 via m3, c3

f4(x)=-m4*x+c4
fit f4(x) '$0.ISF' using 1:6 via m4, c4

f5(x)=-m5*x+c5
fit f5(x) '$0.ISF' using 1:7 via m5, c5

f6(x)=-m6*x+c6
fit f6(x) '$0.ISF' using 1:8 via m6, c6

f7(x)=-m7*x+c7
fit f7(x) '$0.ISF' using 1:9 via m7, c7

f8(x)=-m8*x+c8
fit f8(x) '$0.ISF' using 1:10 via m8, c8

f9(x)=-m9*x+c9
fit f9(x) '$0.ISF' using 1:11 via m9, c9

f10(x)=-m10*x+c10
fit f10(x) '$0.ISF' using 1:12 via m10, c10

f11(x)=-m11*x+c11
fit f11(x) '$0.ISF' using 1:13 via m11, c11

f12(x)=-m12*x+c12
fit f12(x) '$0.ISF' using 1:14 via m12, c12

f13(x)=-m13*x+c13
fit f13(x) '$0.ISF' using 1:15 via m13, c13

f14(x)=-m14*x+c14
fit f14(x) '$0.ISF' using 1:16 via m14, c14

f15(x)=-m15*x+c15
fit f15(x) '$0.ISF' using 1:17 via m15, c15

update '$0.prm'