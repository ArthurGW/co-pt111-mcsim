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

f16(x)=-m16*x+c16
fit f16(x) '$0.ISF' using 1:18 via m16, c16

f17(x)=-m17*x+c17
fit f17(x) '$0.ISF' using 1:19 via m17, c17

f18(x)=-m18*x+c18
fit f18(x) '$0.ISF' using 1:20 via m18, c18

f19(x)=-m19*x+c19
fit f19(x) '$0.ISF' using 1:21 via m19, c19

f20(x)=-m20*x+c20
fit f20(x) '$0.ISF' using 1:22 via m20, c20

f21(x)=-m21*x+c21
fit f21(x) '$0.ISF' using 1:23 via m21, c21

f22(x)=-m22*x+c22
fit f22(x) '$0.ISF' using 1:24 via m22, c22

f23(x)=-m23*x+c23
fit f23(x) '$0.ISF' using 1:25 via m23, c23

f24(x)=-m24*x+c24
fit f24(x) '$0.ISF' using 1:26 via m24, c24

f25(x)=-m25*x+c25
fit f25(x) '$0.ISF' using 1:27 via m25, c25

f26(x)=-m26*x+c26
fit f26(x) '$0.ISF' using 1:28 via m26, c26

f27(x)=-m27*x+c27
fit f27(x) '$0.ISF' using 1:29 via m27, c27

f28(x)=-m28*x+c28
fit f28(x) '$0.ISF' using 1:30 via m28, c28

f29(x)=-m29*x+c29
fit f29(x) '$0.ISF' using 1:31 via m29, c29

f30(x)=-m30*x+c30
fit f30(x) '$0.ISF' using 1:32 via m30, c30

f31(x)=-m31*x+c31
fit f31(x) '$0.ISF' using 1:33 via m31, c31

update '$0.prm'