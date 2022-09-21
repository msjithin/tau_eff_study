

void mutau_analyzer::selections(float weight, int shift, string uncObject)
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
	std::vector<int> event_mu, event_tau;
	event_mu.clear();
	event_tau.clear();
	muCand.clear();
	tauCand.clear();
	aisrtauCand.clear();
	jetCand.clear();
	for (int i = 0; i < nJet; i++)
		jetPt->at(i) = orginal_jetPt[i];

	int mu_index = -1;
	int tau_index = -1;

	if (is_MC)
		event_weight = weight;
	else
		event_weight = 1.0;

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
			gentau1.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
		}
		if(abs(mcMotherPID->at(i))==pdgid_tocheck && mcPID->at(i)==-15  )
		{
			gentau2.SetPtEtaPhiE(mcPt->at(i) , mcEta->at(i) , mcPhi->at(i)  , mcE->at(i));
		}
	}
	//cout<< endl;

	// cout<<"P4 sum  "<< genhiggs.Pt() << " " << genhiggs.Eta() << " " << genhiggs.Phi() <<" " << genhiggs.E() <<endl;
	// cout<< "muon size "<< genmucand.size() <<" tau size "<< gentaucand.size() <<endl;
	// cout<< endl;
	muCand = simple_mu_cand(20, 2.4, 0); ///// ele selected
	tauCand = simple_tau_cand(30.0, 2.3, 0);
	metP4.SetPtEtaPhiE(pfMET, 0, pfMETPhi, pfMET);

	plot_boosted = false;
	event_met_p4 = metP4;
	// muCand.clear();  tauCand.clear();
	// muCand = simple_mu_cand(20,2.4, shift);  ///// muons selected
	mu_index = -1;
	tau_index = -1;
	if (nMu>0 && nTau > 0)
	{

		pair<int, int> selected_indices = get_index();
		make_iso_plot = false;
		mu_index = selected_indices.first;
		tau_index = selected_indices.second;
		if (mu_index >= 0 && tau_index >= 0)
		{
			//setMyEleTau(mu_index, tau_index, metP4, shift);
			MuIndex = mu_index;
			TauIndex = tau_index;
			my_muP4.SetPtEtaPhiE(muPt->at(MuIndex), muEta->at(MuIndex),
								muPhi->at(MuIndex), muE->at(MuIndex));
			my_tauP4.SetPtEtaPhiE(tau_Pt->at(TauIndex), tau_Eta->at(TauIndex), tau_Phi->at(TauIndex), tau_Energy->at(TauIndex));

			my_genmatching_l1 = myGenMaching1(MuIndex);
			my_genmatching_l2 = myGenMaching(TauIndex);
			jetCand.clear();
			jetCand = getJetCand(MuIndex, TauIndex);
			my_njets = jetCand.size();

			applySf = 1.0;
			if (is_MC)
				applySf = getScaleFactors(my_muP4.Pt(),
										  my_tauP4.Pt(),
										  my_muP4.Eta(),
										  my_tauP4.Eta(),
										  tau_DecayMode->at(TauIndex),
										  my_genmatching_l2,
										  false /// this is set to true for fake bakground
				);

			// if(debug)cout<<" sf : "<<getScaleFactors( EleIndex[0] , TauIndex[0] , false , is_MC , debug ) <<endl;
			// cout<<" sf : "<< applySf <<endl;
			event_weight = event_weight * applySf;

			double mvis = (my_muP4 + my_tauP4).M();
			double higgsPt = (my_muP4 + my_tauP4).Pt();
			double tot_tr_mass = (my_muP4 + my_tauP4 + my_metP4).M();
			float deltaR = my_muP4.DeltaR(my_tauP4);

			if (thirdLeptonVeto(MuIndex, TauIndex))
			{

				bool pass_bjet_veto = ((bJet_medium(MuIndex, TauIndex).size() == 0) && (bJet_loose(MuIndex, TauIndex).size() < 2));
				if (pass_bjet_veto)
				{
				  if(my_muP4.DeltaR(my_tauP4) > 0.5)
				    {
					// cout<<__LINE__<<endl;
					plot_resolved_taus(my_muP4, my_tauP4, TauIndex, "6", event_weight);
					if (higgsPt > 65)
					{
						plot_resolved_taus(my_muP4, my_tauP4, TauIndex, "7", event_weight);
						if (mvis < 125)
						{
							plot_resolved_taus(my_muP4, my_tauP4, TauIndex, "8", event_weight);
							if (my_metP4.Pt() > 105)
							{
								// cout<<__LINE__<<endl;
								plot_resolved_taus(my_muP4, my_tauP4, TauIndex, "9", event_weight);
							}
						}
					}
				    }
				}
			}
		}
	}

	muCand.clear();
	tauCand.clear();
	muCand = simple_mu_cand(20, 2.4, 0); ///// ele selected
	tauCand = simple_boostedtau_cand(30.0, 2.3, 0);
	metP4.SetPtEtaPhiE(pfMET, 0, pfMETPhi, pfMET);

	plot_boosted = false;
	event_met_p4 = metP4;
	// muCand.clear();  tauCand.clear();
	// muCand = simple_mu_cand(20,2.4, shift);  ///// muons selected
	mu_index = -1;
	tau_index = -1;

	if (nMu> 0 && nBoostedTau>0 )
	{
		// cout << __LINE__ << endl;
		pair<int, int> selected_indices = get_index_2();
		// cout << __LINE__ << endl;
		make_iso_plot = false;
		mu_index = selected_indices.first;
		tau_index = selected_indices.second;

		// cout << __LINE__ << endl;

		if (mu_index >= 0 && tau_index >= 0)
		{
			// cout << __LINE__ << endl;

			setMyEleTau_boosted(mu_index, tau_index, metP4, shift);
			TauIndex = tau_index;
			MuIndex = mu_index;
			double mvis = (my_muP4 + my_tauP4).M();
			double higgsPt = (my_muP4 + my_tauP4).Pt();
			double tot_tr_mass = (my_muP4 + my_tauP4 + my_metP4).M();
			float deltaR = my_muP4.DeltaR(my_tauP4);

			if (thirdLeptonVeto_boosted(MuIndex, TauIndex))
			{

				bool pass_bjet_veto = ((bJet_medium_boosted(MuIndex, TauIndex).size() == 0) && (bJet_loose_boosted(MuIndex, TauIndex).size() < 2));
				if (pass_bjet_veto)
				{
				  if(deltaR < 0.5)
				    {
					// cout << __LINE__ << endl;
					plot_boosted_taus(my_muP4, my_tauP4, TauIndex, "6", event_weight);
					// cout << __LINE__ << endl;
					if (higgsPt > 65)
					{
						plot_boosted_taus(my_muP4, my_tauP4, TauIndex, "7", event_weight);
						if (mvis < 125)
						{
							plot_boosted_taus(my_muP4, my_tauP4, TauIndex, "8", event_weight);
							if (my_metP4.Pt() > 105)
							{
								// cout<<__LINE__<<endl;
								plot_boosted_taus(my_muP4, my_tauP4, TauIndex, "9", event_weight);
							}
						}
					}
				    }
				}
			}
		}
	}
}

void mutau_analyzer::plot_resolved_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{

	double tauPt = muP4.Pt();
	double subleadingtauPt = tauP4.Pt();
	double deltaR = muP4.DeltaR(tauP4);
	double higgsPt = (muP4 + tauP4).Pt();
	plotFill("tauPt_raw_" + hnumber, tauPt, 970, 30, 1000, event_weight);  // no isolation
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

void mutau_analyzer::plot_boosted_taus(TLorentzVector muP4, TLorentzVector tauP4, int tauindex, string hnumber, double event_weight)
{
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

pair<int, int> mutau_analyzer::get_index()
{

  int mu_index = -1;
  int tau_index = -1;
  float relMuIso, relMuIso_old, relMuIso_v3;

  if (nMu > 0 && nTau>0)
  {
    for (int i = 0; i < nMu; i++)
    {
      for (int j = 0; j < nTau; j++)
      {
        int iMu = i;
        int iTau = j;
        TLorentzVector mu_p4, tau_p4;
        mu_p4.SetPtEtaPhiE(muPt->at(iMu), muEta->at(iMu), muPhi->at(iMu), muE->at(iMu));
        tau_p4.SetPtEtaPhiE(tau_Pt->at(iTau), tau_Eta->at(iTau), tau_Phi->at(iTau), tau_Energy->at(iTau));
        float mu_tau_dr = mu_p4.DeltaR(tau_p4);
        double higgspt = (mu_p4 + tau_p4).Pt();
        relMuIso = (muPFChIso->at(iMu) + max(muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 * muPFPUIso->at(iMu), 0.0) - tau_Pt->at(iTau)) / (muPt->at(iMu));

        bool pass_bjet_veto = ((bJet_medium(iMu, iTau).size() == 0) && (bJet_loose(iMu, iTau).size() < 2));
        bool pass3rdLeptonVeto = (passDiMuonVeto(iMu) == true && eVetoZTTp001dxyz(iMu, iTau) && mVetoZTTp001dxyz(iMu, iTau));
        
        double dr_aa = mu_p4.DeltaR(gentau1);
			  double dr_ab = mu_p4.DeltaR(gentau2);
			  double dr_ba = tau_p4.DeltaR(gentau1);
			  double dr_bb = tau_p4.DeltaR(gentau2);


        if (  muPt->at(iMu) > 20 
			&& fabs(muEta->at(iMu)) < 2.4
			&& fabs(muDz->at(iMu)) < 0.2 
			&& fabs(muD0->at(iMu)) < 0.045
			&& muCharge->at(iMu) * tau_Charge->at(iTau) < 0 
	        && muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
			&& tau_Pt->at(iTau) > 30 
			&& fabs(tau_Eta->at(iTau)) < 2.3
			&& tau_LeadChargedHadron_dz->at(iTau) < 0.2
	        //&& relMuIso < 0.15               // muon relatice isolation
            //&& tau_byMediumDeepTau2017v2p1VSjet->at(iTau)==1
            && tau_IDbits->at(iTau)>>1&1==1
            && ( tau_DecayMode->at(iTau)>=0)
            && ( tau_byVLooseDeepTau2017v2p1VSe->at(iTau)==1 
			&& tau_byTightDeepTau2017v2p1VSmu->at(iTau)==1)
            && (TriggerSelection(mu_p4, tau_p4) == true )
            //  && (thirdLeptonVeto(iMu, iTau))
            //  && (pass_bjet_veto)
            && (mu_tau_dr > 0.5)
            && ( dr_aa < 0.1 ||  dr_ab < 0.1 ||  dr_ba < 0.1 ||  dr_bb < 0.1)
        )
        {
          return make_pair(iMu, iTau);
          /////////// okay we found the mu-au pair, exit these loops
        }
      }
    }
  }

  return make_pair(-1, -1);
}
pair<int, int> mutau_analyzer::get_index_2()
{

  int mu_index = -1;
  int tau_index = -1;
  
  if (nMu > 0 && nTau>0)
  {
    for (int i = 0; i < nMu; i++)
    {
      for (int j = 0; j < nBoostedTau; j++)
      {
        int iMu = i;
        int iTau = j;
        TLorentzVector mu_p4, tau_p4;
        mu_p4.SetPtEtaPhiE(muPt->at(iMu), muEta->at(iMu), muPhi->at(iMu), muE->at(iMu));
        tau_p4.SetPtEtaPhiE(boostedTauPt->at(iTau), boostedTauEta->at(iTau), boostedTauPhi->at(iTau), boostedTauEnergy->at(iTau));
        
        float mu_tau_dr = mu_p4.DeltaR(tau_p4);
        double higgspt = (mu_p4 + tau_p4).Pt();

        double dr_aa = mu_p4.DeltaR(gentau1);
			  double dr_ab = mu_p4.DeltaR(gentau2);
			  double dr_ba = tau_p4.DeltaR(gentau1);
			  double dr_bb = tau_p4.DeltaR(gentau2);

        if ( muPt->at(iMu) > 20 
			&& fabs(muEta->at(iMu)) < 2.4
			&& fabs(muDz->at(iMu)) < 0.2 
			&& fabs(muD0->at(iMu)) < 0.045
			&& muCharge->at(iMu) * boostedTauCharge->at(iTau) < 0 
	        && muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
			&& tau_p4.Pt() > 30
			&& fabs(tau_p4.Eta()) < 2.3
			&& boostedTauZImpact->at(iTau) < 0.2
            && muIDbit->at(iMu) >> 1 & 1 == 1 // muon id
            && boostedTaupfTausDiscriminationByDecayModeFindingNewDMs->at(iTau) > 0.5
            && boostedTauByMVA6LooseElectronRejection->at(iTau) > 0.5
            && boostedTauByTightMuonRejection3->at(iTau) > 0.5
            && mu_tau_dr < 0.5 
            && ( dr_aa < 0.1 ||  dr_ab < 0.1 ||  dr_ba < 0.1 ||  dr_bb < 0.1)
        )
        {
          return make_pair(iMu, iTau);
        }
      }
    }
  }

  return make_pair(-1, -1);
}