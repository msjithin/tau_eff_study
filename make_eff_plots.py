import ROOT as R
from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin
from array import array
from itertools import product

R.gStyle.SetFrameLineWidth(1)
R.gStyle.SetLineWidth(2)
R.gStyle.SetOptStat(0)
R.gStyle.SetEndErrorSize(5)

def make_legend():
  output = R.TLegend(0.85, 0.2, 0.99, 0.50, "", "brNDC")
  #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
  output.SetLineWidth(1)
  output.SetLineStyle(1)
  output.SetFillStyle(1001) #0
  output.SetFillColor(0)
  output.SetBorderSize(1)
  output.SetTextFont(42)
  output.SetNColumns(1)                                                                                                                                                                                                                                              
                                       
  return output

def add_lumi():
	lowX=0.550
	lowY=0.82
	lumi  = R.TPaveText(lowX, lowY+0.06, lowX+0.40, lowY+0.15, "NDC")
	lumi.SetBorderSize(   0 )
	lumi.SetFillStyle(    0 )
	lumi.SetTextAlign(   32 )#12
	lumi.SetTextColor(    1 )
	lumi.SetTextSize(0.03)
	lumi.SetTextFont (   42 )
	lumiProcessed="41.52"
	channel_ = 'tautau'
	if channel_=="combined":
		lumi.AddText("4 channels combined 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="mutau":
		lumi.AddText("#mu#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="etau":
		lumi.AddText("e#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="tautau":
		lumi.AddText("#tau_{h}#tau_{h} 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	if channel_=="emu":
		lumi.AddText("e#mu 2017, "+lumiProcessed+" fb^{-1} (13 TeV)")
	return lumi

mzp_map = {'2' : 'MZp_1000_MChi_1' , '7': 'MZp_100_MChi_1', '10': 'MZp_1500_MChi_1',  'DY': 'DYJetsToLL'}
canv=R.TCanvas("canvas","",0,0,1300,1200)
canv.cd()

vars = ['HiggsPt', 'tauPt', 'deltaR', 'subleadingtauPt']
selections = ['6', '7', '8', '9']
idtypes = { 'deep' : ['VVVLoose', 'Loose', 'Medium', 'Tight', 'VVTight'] ,
            'boosted' : ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight'] ,
            '2017'    : ['VLoose' , 'Loose', 'Medium', 'Tight', 'VTight'] ,
            '2016'    : ['VLoose' , 'Loose', 'Medium', 'Tight', 'VTight']
   }
new_binning = array('d', [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 300, 350, 400, 450,  500, 600, 700,800,900, 1000])


parameters = list(product(vars, selections, idtypes))

for idx in mzp_map:
	mzp_idx = idx
	file_name = 'output/Zpbaryonic2017_'+mzp_idx+'.root'
	if idx=='DY':
		file_name = 'output/DYJetsToLL_M-50_TuneCP5.root'
	infile = R.TFile(file_name, 'update')

	for var, selection, idtype in parameters:
		signal = mzp_map[mzp_idx]

		raw_name = ''
		
		if idtype == 'boosted':
			raw_name = var+'_boostedraw_'+selection
		else:
			raw_name = var+'_raw_'+selection

		# print raw_name
		hist_raw = infile.Get(raw_name)
		idlist = idtypes[idtype]

		hist_1 = infile.Get(var+'_'+idtype+idlist[0]+'_'+selection)
		hist_2 = infile.Get(var+'_'+idtype+idlist[1]+'_'+selection)
		hist_3 = infile.Get(var+'_'+idtype+idlist[2]+'_'+selection)
		hist_4 = infile.Get(var+'_'+idtype+idlist[3]+'_'+selection)
		hist_5 = infile.Get(var+'_'+idtype+idlist[4]+'_'+selection)


		if var != 'deltaR':
			tmpHist_num1 = hist_1.Rebin(22, var+'_'+idtype+idlist[0]+'_'+selection, new_binning)
			tmpHist_num2 = hist_2.Rebin(22, var+'_'+idtype+idlist[1]+'_'+selection, new_binning)
			tmpHist_num3 = hist_3.Rebin(22, var+'_'+idtype+idlist[2]+'_'+selection, new_binning)
			tmpHist_num4 = hist_4.Rebin(22, var+'_'+idtype+idlist[3]+'_'+selection, new_binning)
			tmpHist_num5 = hist_5.Rebin(22, var+'_'+idtype+idlist[4]+'_'+selection, new_binning)
			tmpHist_den = hist_raw.Rebin(22, raw_name, new_binning )  
		else:
			dr_bins = 5
			tmpHist_num1 = hist_1.Rebin(dr_bins, var+'_'+idtype+idlist[0]+'_'+selection)
			tmpHist_num2 = hist_2.Rebin(dr_bins, var+'_'+idtype+idlist[1]+'_'+selection)
			tmpHist_num3 = hist_3.Rebin(dr_bins, var+'_'+idtype+idlist[2]+'_'+selection)
			tmpHist_num4 = hist_4.Rebin(dr_bins, var+'_'+idtype+idlist[3]+'_'+selection)
			tmpHist_num5 = hist_5.Rebin(dr_bins, var+'_'+idtype+idlist[4]+'_'+selection)
			tmpHist_den = hist_raw.Rebin(dr_bins, raw_name)  		

		tmpHist_num1.Divide(tmpHist_den)
		tmpHist_num2.Divide(tmpHist_den)
		tmpHist_num3.Divide(tmpHist_den)
		tmpHist_num4.Divide(tmpHist_den)
		tmpHist_num5.Divide(tmpHist_den)

		if var=='deltaR' and idtype=='boosted':		
			tmpHist_num1.GetXaxis().SetRangeUser(0, 1.0) 
			tmpHist_num2.GetXaxis().SetRangeUser(0, 1.0) 
			tmpHist_num3.GetXaxis().SetRangeUser(0, 1.0) 
			tmpHist_num4.GetXaxis().SetRangeUser(0, 1.0) 
			tmpHist_num5.GetXaxis().SetRangeUser(0, 1.0) 
		elif var=='deltaR' and idtype!='boosted':		
			tmpHist_num1.GetXaxis().SetRangeUser(0, 4.0) 
			tmpHist_num2.GetXaxis().SetRangeUser(0, 4.0) 
			tmpHist_num3.GetXaxis().SetRangeUser(0, 4.0) 
			tmpHist_num4.GetXaxis().SetRangeUser(0, 4.0) 
			tmpHist_num5.GetXaxis().SetRangeUser(0, 4.0) 

		tmpHist_num1.GetXaxis().SetTitle(var)

		# if var=='tauPt':
		# 	var = 'muonPt'
		# elif var=='subleadingtauPt':
		# 	var = 'tauPt'
		tmpHist_num1.SetTitle(var+ ' efficiency '+idtype+' tau isolation '+signal)
		tmpHist_num1.SetMaximum(1.2)
		tmpHist_num1.SetMinimum(0.0)

		tmpHist_num1.SetMarkerStyle(20)
		tmpHist_num2.SetMarkerStyle(21)
		tmpHist_num3.SetMarkerStyle(22)
		tmpHist_num4.SetMarkerStyle(23)
		tmpHist_num5.SetMarkerStyle(33)

		tmpHist_num1.SetMarkerSize(3)
		tmpHist_num2.SetMarkerSize(3)
		tmpHist_num3.SetMarkerSize(3)
		tmpHist_num4.SetMarkerSize(3)
		tmpHist_num5.SetMarkerSize(3)

		tmpHist_num1.SetMarkerColor(2)
		tmpHist_num2.SetMarkerColor(3)
		tmpHist_num3.SetMarkerColor(4)
		tmpHist_num4.SetMarkerColor(6)
		tmpHist_num5.SetMarkerColor(7)

		tmpHist_num1.SetLineColor(2)
		tmpHist_num2.SetLineColor(3)
		tmpHist_num3.SetLineColor(4)
		tmpHist_num4.SetLineColor(6)
		tmpHist_num5.SetLineColor(7)

		nbinsX = tmpHist_num1.GetNbinsX()

		tmpHist_num1.Draw("e1")
		tmpHist_num2.Draw("e1same")
		tmpHist_num3.Draw("e1same")
		tmpHist_num4.Draw("e1same")
		tmpHist_num5.Draw("e1same")

		legende=make_legend()
		legende.AddEntry(tmpHist_num1,idlist[0],"lp")
		legende.AddEntry(tmpHist_num2,idlist[1],"lp")
		legende.AddEntry(tmpHist_num3,idlist[2],"lp")
		legende.AddEntry(tmpHist_num4,idlist[3],"lp")
		legende.AddEntry(tmpHist_num5,idlist[4],"lp")

		legende.Draw("same")
		l1=add_lumi()
		l1.Draw("same")

		canv.SaveAs('plots/eff_2017tautau_'+idtype+'_'+signal+'_'+var+'_'+selection+'.png')

