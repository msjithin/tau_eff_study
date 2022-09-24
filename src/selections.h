void tautau_analyzer::selections(float weight, int shift, string uncObject)
{
	check_unc = false; // set true for printing unc pt, values

	double event_weight = weight;
	TLorentzVector metP4, event_met_p4;
	if (shift > 0)
		unc_shift = "up";
	else if (shift < 0)
		unc_shift = "down";
	else
		unc_shift = "nominal";
	shift_index = shift;
	selected_systematic = uncObject;
	// cout<<" selected_systematic = "<< selected_systematic << " shift = "<< shift<<endl;
	int event_tau1, event_tau2;
	int index_tau1, index_tau2;
	std::vector<int> tauCand;
	index_tau1 = index_tau1 = -1;
	tau1Cand.clear();
	tau2Cand.clear();
	jetCand.clear();
	for (int i = 0; i < nJet; i++)
		jetPt->at(i) = orginal_jetPt[i];

	vector<int> gentaucand, genmucand, genhiggscand;
	gentaucand.clear();
	genmucand.clear();
	genhiggscand.clear();

	int pdgid_tocheck = 25;
	if(found_DYjet_sample)
		pdgid_tocheck =23;
	for(int i=0; i<nMC; i++){
		// cout<< i << " " << mcPID->at(i) << " " << " " << mcMotherPID->at(i) <<endl;
		if(abs(mcPID->at(i))==pdgid_tocheck ){
			genhiggscand.push_back(i);
		}
	}

	for(int i=0; i<nMC; i++){
		if(abs(mcMotherPID->at(i))==pdgid_tocheck && mcPID->at(i)==15  )
		{
			//cout<<"Found +taus"<<endl;
			genmucand.push_back(i);
			gentau1.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
		}
		if(abs(mcMotherPID->at(i))==pdgid_tocheck && mcPID->at(i)==-15  )
		{
			//cout<<"Found -taus"<<endl;
			gentaucand.push_back(i);
			gentau2.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
		}
	}

	double gentauPt = max( gentau1.Pt(), gentau2.Pt() );
	double gensubleadingtauPt = min( gentau1.Pt(), gentau2.Pt() );
	double gentau_deltaR = gentau1.DeltaR(gentau2);
	double gentau_higgsPt = (gentau1 + gentau2).Pt();
	if( gentau1.Pt()>40 && gentau2.Pt()>40  && abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 )
	{
		plotFill("gentauPt_raw_0", gentauPt, 970, 30, 1000, event_weight);
		plotFill("gensubleadingtauPt_raw_0", gensubleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("genHiggsPt_raw_0", gentau_higgsPt, 970, 30, 1000, event_weight);
		plotFill("gendeltaR_raw_0", gentau_deltaR, 600, 0, 6, event_weight);
	}
	// Setting nominal values
	// index_tau1 = getTau1Cand(40.0, 2.1, 0);
	// index_tau2 = getTau2Cand(40.0, 2.1, index_tau1, 0);
	metP4.SetPtEtaPhiE(pfMET, 0, pfMETPhi, pfMET);

	// tauCand.clear();
	// tauCand = simple_tau_cand(40.0, 2.1, 0);

	if (nTau > 1)
	{
		/// add hist
		// if (shift == 0 && tauCand.size() > 1)
		// 	fillHist("2", tauCand[0], tauCand[1], false, event_weight);
		nevents3++;
		pair<int, int> selected_indices = get_index();
		index_tau1 = selected_indices.first;
		index_tau2 = selected_indices.second;
		if (index_tau1 >= 0 && index_tau2 >= 0)
		{
			if( gentau1.Pt()>40 && gentau2.Pt()>40  && abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 )
			{
				plotFill("gentauPt_raw_1", gentauPt, 970, 30, 1000, event_weight);
				plotFill("gensubleadingtauPt_raw_1", gensubleadingtauPt, 970, 30, 1000, event_weight);
				plotFill("genHiggsPt_raw_1", gentau_higgsPt, 970, 30, 1000, event_weight);
				plotFill("gendeltaR_raw_1", gentau_deltaR, 600, 0, 6, event_weight);
			}
			nevents4++;
			nGoodTauPassed += event_weight;
			//setMyEleTau(index_tau1, index_tau2, metP4, shift);
			Tau1Index = index_tau1;
			Tau2Index = index_tau2;
			my_tau1P4.SetPtEtaPhiE(tau_Pt->at(Tau1Index),tau_Eta->at(Tau1Index) ,tau_Phi->at(Tau1Index), tau_Energy->at(Tau1Index));
  			my_tau2P4.SetPtEtaPhiE(tau_Pt->at(Tau2Index),tau_Eta->at(Tau2Index) ,tau_Phi->at(Tau2Index), tau_Energy->at(Tau2Index));
			my_metP4.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET); 
			// third lepton veto
			if (thirdLeptonVeto(Tau1Index, Tau2Index))
			{
				plot_resolved_taus(my_tau1P4, my_tau2P4, Tau1Index, "5", event_weight);
				// if( my_tau1P4.DeltaR(my_tau2P4) > 0.5)
				{
					plot_resolved_taus(my_tau1P4, my_tau2P4, Tau1Index, "6", event_weight);

					// makeTestPlot("i", 0,0,0,event_weight);
					double mvis = (my_tau1P4 + my_tau2P4).M();
					double higgsPt = (my_tau1P4 + my_tau2P4).Pt();
					if (higgsPt > 65)
					{
						nevents6++;
						plot_resolved_taus(my_tau1P4, my_tau2P4, Tau1Index, "7", event_weight);

						if (mvis < 125)
						{
							plot_resolved_taus(my_tau1P4, my_tau2P4, Tau1Index, "8", event_weight);

							if (my_metP4.Pt() > 105)
							{
								plot_resolved_taus(my_tau1P4, my_tau2P4, Tau1Index, "9", event_weight);
							}
						}
					}
				}
			}
		}
	}
	
	
	//// cout << __LINE__ << endl;
	// tauCand.clear();
	// tauCand = simple_tau_cand_2(40.0, 2.1, 0);
	if (nBoostedTau> 1)
	{
		//// cout << __LINE__ << endl;
		/// add hist
		// if (shift == 0 && tauCand.size() > 1)
		// 	fillHist("2", tauCand[0], tauCand[1], false, event_weight);
		nevents3++;
		pair<int, int> selected_indices = get_index_2();
		index_tau1 = selected_indices.first;
		index_tau2 = selected_indices.second;
		if (index_tau1 >= 0 && index_tau2 >= 0)
		{
			// cout << __LINE__ << endl;

			nevents4++;
			nGoodTauPassed += event_weight;
			if( gentau1.Pt()>40 && gentau2.Pt()>40  && abs(gentau1.Eta())<2.3 && abs(gentau2.Eta())<2.3 )
			{
				plotFill("gentauPt_boostedraw_1", gentauPt, 970, 30, 1000, event_weight);
				plotFill("gensubleadingtauPt_boostedraw_1", gensubleadingtauPt, 970, 30, 1000, event_weight);
				plotFill("genHiggsPt_boostedraw_1", gentau_higgsPt, 970, 30, 1000, event_weight);
				plotFill("gendeltaR_boostedraw_1", gentau_deltaR, 600, 0, 6, event_weight);
			}
			// setMyEleTau(index_tau1, index_tau2, metP4, shift);
			Tau1Index = index_tau1;
			Tau2Index = index_tau2;
			// // cout << __LINE__ << endl;

			my_tau1P4.SetPtEtaPhiE(boostedTauPt->at(index_tau1), boostedTauEta->at(index_tau1), boostedTauPhi->at(index_tau1), boostedTauEnergy->at(index_tau1));
			my_tau2P4.SetPtEtaPhiE(boostedTauPt->at(index_tau2), boostedTauEta->at(index_tau2), boostedTauPhi->at(index_tau2), boostedTauEnergy->at(index_tau2));
			my_metP4.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET); 
			if (thirdLeptonVeto_boosted(Tau1Index, Tau2Index))
			{
				// cout << __LINE__ << endl;
				plot_boosted_taus(my_tau1P4, my_tau2P4, Tau1Index, "5", event_weight);
				// if( my_tau1P4.DeltaR(my_tau2P4) < 0.5)
				{
					plot_boosted_taus(my_tau1P4, my_tau2P4, Tau1Index, "6", event_weight);

					// makeTestPlot("i", 0,0,0,event_weight);
					double mvis = (my_tau1P4 + my_tau2P4).M();
					double higgsPt = (my_tau1P4 + my_tau2P4).Pt();
					if (higgsPt > 65)
					{
						// cout << __LINE__ << endl;

						nevents6++;
						plot_boosted_taus(my_tau1P4, my_tau2P4, Tau1Index, "7", event_weight);

						if (mvis < 125)
						{
							// cout << __LINE__ << endl;

							plot_boosted_taus(my_tau1P4, my_tau2P4, Tau1Index, "8", event_weight);

							if (my_metP4.Pt() > 105)
							{
								plot_boosted_taus(my_tau1P4, my_tau2P4, Tau1Index, "9", event_weight);
							}
						}
					}
				}
			}
		}
	}
}

void tautau_analyzer::plot_resolved_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{

	double gentauPt = max( gentau1.Pt(), gentau2.Pt() );
	double gensubleadingtauPt = min( gentau1.Pt(), gentau2.Pt() );
	double gentau_deltaR = gentau1.DeltaR(gentau2);
	double gentau_higgsPt = (gentau1 + gentau2).Pt();
	plotFill("gentauPt_raw_"+hnumber, gentauPt, 970, 30, 1000, event_weight);
	plotFill("gensubleadingtauPt_raw_"+hnumber, gensubleadingtauPt, 970, 30, 1000, event_weight);
	plotFill("genHiggsPt_raw_"+hnumber, gentau_higgsPt, 970, 30, 1000, event_weight);
	plotFill("gendeltaR_raw_"+hnumber, gentau_deltaR, 600, 0, 6, event_weight);

	double tauPt = muP4.Pt();
	double subleadingtauPt = tauP4.Pt();
	double deltaR = muP4.DeltaR(tauP4);
	double higgsPt = (muP4 + tauP4).Pt();
	plotFill("tauPt_raw_" + hnumber, tauPt, 970, 30, 1000, event_weight);
	plotFill("subleadingtauPt_raw_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
	plotFill("HiggsPt_raw_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	plotFill("deltaR_raw_" + hnumber, deltaR, 600, 0, 6, event_weight);
	if (tau_byVVVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepVVVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepVVVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepVVVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepVVVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (tau_byVVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepVVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepVVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepVVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepVVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (tau_byVLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_byLooseDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_byMediumDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepMedium_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepMedium_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepMedium_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepMedium_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_byTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_byVTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepVTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepVTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepVTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepVTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_byVVTightDeepTau2017v2p1VSjet->at(tauindex) == 1)
	{
		plotFill("tauPt_deepVVTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_deepVVTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_deepVVTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_deepVVTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}

	if (tau_IDbits->at(tauindex) >> 12 & 1 == 1)
	{
		plotFill("tauPt_2017VVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017VVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017VVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017VVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 13 & 1 == 1)
	{
		plotFill("tauPt_2017VLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017VLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017VLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017VLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 14 & 1 == 1)
	{
		plotFill("tauPt_2017Loose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017Loose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017Loose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017Loose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 15 & 1 == 1)
	{
		plotFill("tauPt_2017Medium_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017Medium_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017Medium_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017Medium_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 16 & 1 == 1)
	{
		plotFill("tauPt_2017Tight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017Tight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017Tight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017Tight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 17 & 1 == 1)
	{
		plotFill("tauPt_2017VTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017VTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017VTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017VTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 18 & 1 == 1)
	{
		plotFill("tauPt_2017VVTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2017VVTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2017VVTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2017VVTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 20 & 1 == 1)
	{
		plotFill("tauPt_2016VVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016VVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016VVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016VVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (tau_IDbits->at(tauindex) >> 21 & 1 == 1)
	{
		plotFill("tauPt_2016VLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016VLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016VLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016VLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}

	if (tau_IDbits->at(tauindex) >> 22 & 1 == 1)
	{
		plotFill("tauPt_2016Loose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016Loose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016Loose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016Loose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 23 & 1 == 1)
	{
		plotFill("tauPt_2016Medium_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016Medium_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016Medium_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016Medium_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (tau_IDbits->at(tauindex) >> 24 & 1 == 1)
	{
		plotFill("tauPt_2016Tight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016Tight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016Tight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016Tight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (tau_IDbits->at(tauindex) >> 25 & 1 == 1)
	{
		plotFill("tauPt_2016VTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016VTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016VTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016VTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	if (tau_IDbits->at(tauindex) >> 26 & 1 == 1)
	{
		plotFill("tauPt_2016VVTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_2016VVTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_2016VVTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_2016VVTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
}

void tautau_analyzer::plot_boosted_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{

	double gentauPt = max( gentau1.Pt(), gentau2.Pt() );
	double gensubleadingtauPt = min( gentau1.Pt(), gentau2.Pt() );
	double gentau_deltaR = gentau1.DeltaR(gentau2);
	double gentau_higgsPt = (gentau1 + gentau2).Pt();
	plotFill("gentauPt_raw_"+hnumber, gentauPt, 970, 30, 1000, event_weight);
	plotFill("gensubleadingtauPt_raw_"+hnumber, gensubleadingtauPt, 970, 30, 1000, event_weight);
	plotFill("genHiggsPt_raw_"+hnumber, gentau_higgsPt, 970, 30, 1000, event_weight);
	plotFill("gendeltaR_raw_"+hnumber, gentau_deltaR, 600, 0, 6, event_weight);

	double tauPt = muP4.Pt();
	double subleadingtauPt = tauP4.Pt();
	double deltaR = muP4.DeltaR(tauP4);
	double higgsPt = (muP4 + tauP4).Pt();

	// cout<<__LINE__<<endl;
	plotFill("tauPt_boostedraw_" + hnumber, tauPt, 970, 30, 1000, event_weight);
	plotFill("subleadingtauPt_boostedraw_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
	plotFill("deltaR_boostedraw_" + hnumber, deltaR, 600, 0, 6, event_weight);
	plotFill("HiggsPt_boostedraw_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	// cout<<__LINE__<<endl;
	if (boostedTauByVLooseIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
	{
		plotFill("tauPt_boostedVLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_boostedVLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_boostedVLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_boostedVLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);
	}
	// cout<<__LINE__<<endl;
	if (boostedTauByLooseIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
	{
		plotFill("tauPt_boostedLoose_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_boostedLoose_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_boostedLoose_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_boostedLoose_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	// cout<<__LINE__<<endl;
	if (boostedTauByMediumIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
	{
		plotFill("tauPt_boostedMedium_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_boostedMedium_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_boostedMedium_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_boostedMedium_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (boostedTauByTightIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
	{
		plotFill("tauPt_boostedTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_boostedTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_boostedTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_boostedTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
	if (boostedTauByVTightIsolationMVArun2v1DBoldDMwLTNew->at(tauindex) == 1)
	{
		plotFill("tauPt_boostedVTight_" + hnumber, tauPt, 970, 30, 1000, event_weight);
		plotFill("subleadingtauPt_boostedVTight_" + hnumber, subleadingtauPt, 970, 30, 1000, event_weight);
		plotFill("deltaR_boostedVTight_" + hnumber, deltaR, 600, 0, 6, event_weight);
		plotFill("HiggsPt_boostedVTight_" + hnumber, higgsPt, 970, 30, 1000, event_weight);

	}
}


pair<int, int> tautau_analyzer::get_index()
{

  int mu_index = -1;
  int tau_index = -1;
  float relMuIso, relMuIso_old, relMuIso_v3;
  bool plot_boosted = false;
  bool is_boosted = false;


  if ( nTau > 1)
  {
    for (int i = 0; i < nTau; i++)
    {
      for (int j = i+1; j < nTau; j++)
      {
        int iTau1 = i;
        int iTau2 = j;
        TLorentzVector tau1_p4, tau2_p4;
        tau1_p4.SetPtEtaPhiE(tau_Pt->at(iTau1), tau_Eta->at(iTau1), tau_Phi->at(iTau1), tau_Energy->at(iTau1));
        tau2_p4.SetPtEtaPhiE(tau_Pt->at(iTau2), tau_Eta->at(iTau2), tau_Phi->at(iTau2), tau_Energy->at(iTau2));
        float tau_tau_dr = tau1_p4.DeltaR(tau2_p4);
        double higgspt = (tau1_p4 + tau2_p4).Pt();
        
        double dr_aa = tau1_p4.DeltaR(gentau1);
			  double dr_ab = tau1_p4.DeltaR(gentau2);
			  double dr_ba = tau2_p4.DeltaR(gentau1);
			  double dr_bb = tau2_p4.DeltaR(gentau2);

        bool pass3rdLeptonVeto = ( eVetoZTTp001dxyz(iTau1, iTau2) && mVetoZTTp001dxyz(iTau1, iTau2) );
        if (   tau1_p4.Pt() > 40 && tau2_p4.Pt() > 40
				&& abs(tau1_p4.Eta()) < 2.3 && abs(tau2_p4.Eta()) < 2.3
				&& tau_Charge->at(iTau1) * tau_Charge->at(iTau2) < 0 
            // && tau_byMediumDeepTau2017v2p1VSjet->at(iTau1) == 1
            // && tau_byMediumDeepTau2017v2p1VSjet->at(iTau2) == 1
            && tau_DecayMode->at(iTau1) >= 0
            && tau_DecayMode->at(iTau2) >= 0
            && (tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau1) == 1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau1) == 1)
            && (tau_byVVVLooseDeepTau2017v2p1VSe->at(iTau2) == 1 && tau_byVLooseDeepTau2017v2p1VSmu->at(iTau2) == 1)
            && (tau_IDbits->at(iTau1) >> 1 & 1 == 1)
            && (tau_IDbits->at(iTau2) >> 1 & 1 == 1)
            // && (TriggerSelection(tau1_p4, tau2_p4) == true )
		        // && (thirdLeptonVeto(iTau1, iTau2))
            && ( dr_aa < 0.1 ||  dr_ab < 0.1 ||  dr_ba < 0.1 ||  dr_bb < 0.1)
            // && (tau_tau_dr > 0.5)
        )
        {
          return make_pair(iTau1, iTau2);
        }
      }
    }
  }

  return make_pair(-1, -1);
}

pair<int, int> tautau_analyzer::get_index_2()
{

  int mu_index = -1;
  int tau_index = -1;
  float relMuIso, relMuIso_old, relMuIso_v3;
  bool plot_boosted = false;
  bool is_boosted = false;


  if ( nBoostedTau > 1)
  {
    for (int i = 0; i < nBoostedTau; i++)
    {
      for (int j = i+1; j < nBoostedTau; j++)
      {
        int iTau1 = i;
        int iTau2 = j; 
        TLorentzVector tau1_p4, tau2_p4;
        tau1_p4.SetPtEtaPhiE(boostedTauPt->at(iTau1), boostedTauEta->at(iTau1), boostedTauPhi->at(iTau1), boostedTauEnergy->at(iTau1));
        tau2_p4.SetPtEtaPhiE(boostedTauPt->at(iTau2), boostedTauEta->at(iTau2), boostedTauPhi->at(iTau2), boostedTauEnergy->at(iTau2));
        float tau_tau_dr = tau1_p4.DeltaR(tau2_p4);
        double higgspt = (tau1_p4 + tau2_p4).Pt();

        double dr_aa = tau1_p4.DeltaR(gentau1);
			  double dr_ab = tau1_p4.DeltaR(gentau2);
			  double dr_ba = tau2_p4.DeltaR(gentau1);
			  double dr_bb = tau2_p4.DeltaR(gentau2);

        //bool pass3rdLeptonVeto = ( eVetoZTTp001dxyz(iTau1, iTau2) && mVetoZTTp001dxyz(iTau1, iTau2) );
        if ( tau1_p4.Pt() > 40 && tau2_p4.Pt() > 40
				&& abs(tau1_p4.Eta()) < 2.3 && abs(tau2_p4.Eta()) < 2.3  
				&& boostedTauCharge->at(iTau1) * boostedTauCharge->at(iTau2) < 0 
              && boostedTaupfTausDiscriminationByDecayModeFindingNewDMs->at(iTau1) > 0.5
              && boostedTaupfTausDiscriminationByDecayModeFindingNewDMs->at(iTau2) > 0.5
              && boostedTauByLooseMuonRejection3->at(iTau1) > 0.5
              && boostedTauByLooseMuonRejection3->at(iTau2) > 0.5
              && boostedTauDecayMode->at(iTau1) >= 0
              && boostedTauDecayMode->at(iTau2) >= 0
              // && tau_tau_dr < 0.5 
              && ( dr_aa < 0.1 ||  dr_ab < 0.1 ||  dr_ba < 0.1 ||  dr_bb < 0.1)
            )
        {
          return make_pair(iTau1, iTau2);
        }
      }
    }
  }

  return make_pair(-1, -1);
}