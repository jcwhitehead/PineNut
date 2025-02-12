# BEGIN PLOT /MC_TTBAR_PN/.*
LegendAlign=r
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*N
Title=Number of generated events
XLabel=
YLabel=$N$
LogY=0
XmajorCustomMajorTicks=0.5	$\quad$
#XMajorTickMarks=20
XMinorTickMarks=0
LegendXPos=0.05
LegendYPos=0.15
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*pmXS
Title=Fraction of positive and negative weighted events
XLabel=$\mathrm{sgn}(w)$
XCustomMajorTicks=-0.5	$w<0$	0.5	$w\geq0$
YLabel=$\mathrm{d}\sigma/\mathrm{d~sgn}(w)$ [pb]
LogY=0
ShowZero=0
XMinorTickMarks=0
XMin=-1
XMax=1
LegendXPos=0.05
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_XS
Title=Cross-section
XLabel=$w$
YLabel=$\sigma_{\mathrm{tot}}$ [pb]
LogY=0
# END PLOT


# BEGIN PLOT /MC_TTBAR_PN/.*pmN
Title=Number of positive and negative weighted events
XLabel=$\mathrm{sgn}(w)$
YLabel=$N_\pm$
LogY=0
ShowZero=0
XCustomMajorTicks=-0.5	$w<0$	0.5	$w\geq0$
XMinorTickMarks=0
XMin=-1
XMax=1
LegendXPos=0.05
# END PLOT


# BEGIN PLOT /MC_TTBAR_PN/.*tt_absrap
XLabel=$|y^{t\bar{t}}|$
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*tt_m
XLabel=$|m^{t\bar{t}}|$
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*tt_m_vs_tt_absrap
XLabel=$|m^{t\bar{t}}| vs |eta| 0.0, 0.3, 0.6, 0.9, 1.3, 2.5 $
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*tt_pt
XLabel=$|p_{T}^{t\bar{t}}|$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*t_absrap
XLabel=$|y^{t}|$
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_t_pt
XLabel=$p_{T}^{t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_t_pt_vs_t_absrap
XLabel=$p_{T}^{t}$ [GeV] vs |t_eta| 0, 0.4, 0.8, 1.2, 1.6, 2.5
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_top_pt
XLabel=$p_{T}^{t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_top_rap
XLabel=$y^{t}$
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_top_phi
XLabel=$\phi^{t}$
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_top_e
XLabel=$E^{t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_antitop_pt
XLabel=$p_{T}^{\bar{t}}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_antitop_rap
XLabel=$y^{\bar{t}}$
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_antitop_phi
XLabel=$\phi^{\bar{t}}$
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_antitop_e
XLabel=$E^{t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*tt_dR
XLabel=$\Delta R(t,\bar{t})$
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_leadtop_pt
XLabel=$p_{T}^{lead. t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_leadtop_absrap
XLabel=$|y^{lead. t}|$
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_leadtop_phi
XLabel=$\phi^{lead. t}$
YLabel=$\mathrm{Events / Norm.}$
LogY=1
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT

# BEGIN PLOT /MC_TTBAR_PN/.*_leadtop_e
XLabel=$E^{lead. t}$ [GeV]
YLabel=$\mathrm{Events / Norm.}$
LogY=0
RatioPlotYMin=0.9
RatioPlotYMax=1.1
# RatioPlotYCustomMajorTicks=0.999	0.999	1	1	1.001	1.001
# RatioPlotYCustomMinorTicks=0.9992	0.9994	0.9996	0.9998	1.0002	1.0004	1.0006	1.0008
# RatioPlotYMin=0.999
# RatioPlotYMax=1.001
# RatioPlotYMajorTickMarks=10
# RatioPlotYMinorTickMarks=10
# END PLOT
