//
//  BeerLambertUtil.swift
//  Blueberry
//
//  Created by John David Chibuk on 2021-12-12.
//  Copyright © 2021 Blueberry. All rights reserved.
//

import Foundation
import Accelerate

public struct CommonBeerLambert {
    
    //max double 1.7976931348623157E308 / 8192000 (max Int for diffusion) = 2.19444962752e+301 samples
    //
    private static var sumTotalHbOShort = 0.0;
    private static var sumTotalHbRShort = 0.0;
    private static var countTotalShort = 0.0;
    private static var sumTotalHbOLong = 0.0;
    private static var sumTotalHbRLong = 0.0;
    private static var countTotalLong = 0.0;
    
    private static var L1 = 10.3 // mm
    private static var L2 = 27.0 // mm

    //A = L * absorptionCoeff

    //extinction coefficients
    //https://www.researchgate.net/publication/7910860_Extinction_coefficients_of_hemoglobin_for_near-infrared_spectroscopy_of_tissue
    //Wavelength 1 = 740nm
    private static var alphaHbR1 = 1.5
    private static var alphaHbO1 = 0.45
    //Wavelength 2 = 850nm
    private static var alphaHbR2 = 1.1
    private static var alphaHbO2 = 0.75

    //740nm HbR = 1.5    //850nm HbR = 1.1f    //660nm HbR = 1.7f    //880nm HbR = 0.8f     //940nm HbR = 0.9f
    //740nm HbO = 0.45   //850nm HbO = 0.75f   //660nm HbO = 0.35f   //880nm HbO = 1.25f    //940nm HbO = 1.35f

    //DPF value
    //https://www.nature.com/articles/pr19962544/tables/3
    private static var age = 31.0

    //volatile float DPF1 = 4.67f;    //880nm
    //volatile float DPF2 = 4.2f;    //940nm
    //volatile float DPF3 = 6.0f;    //~660nm
    //5.3 + pow(0.049 * age, 0.723); //~720nm

    private static var DPF1 = 6.65    //740nm
    private static var DPF2 = 6.12    //850nm
    
    // 6.650 = 740nm
    // 6.208 = 940nm
    // 6.120 = 850nm
    
    // 4.5 + pow((0.058*age),0.823);    //850nm
    // 4.67 + pow(0.062*age, 0.819);    //~880nm
    // 4.2 + pow((0.07*age),0.9);       //940nm
    // 6.0 + pow(0.03 * age, 0.65);     //~660nm
    // 5.3 + pow(0.049 * age, 0.723); //~720nm

    //    deltaOD1 = log(Ib / It)
    //    deltaOD1 = log(longIR / shortIR)
    //    Ib = baseline light
    //    It = transient light

    public static func changeHbO(_ PD1_IR1: UInt32, PD2_IR1: UInt32, PD1_IR2: UInt32, PD2_IR2:UInt32) -> Double{

        let deltaOD1 = log2(Double(PD2_IR1) / Double(PD1_IR1));
        let deltaOD2 = log2(Double(PD2_IR2) / Double(PD1_IR2));

        let changeHbOVal = (-1)*(alphaHbR1 * (deltaOD2 / Double(DPF2)) - alphaHbO1 * (deltaOD1 / Double(DPF1))) / (L2*(alphaHbR2 * alphaHbO1 - alphaHbR2 * alphaHbO2));
        
        return changeHbOVal;
        
    }

    public static func changeHbR(_ PD1_IR1: UInt32, PD2_IR1: UInt32, PD1_IR2: UInt32, PD2_IR2:UInt32) -> Double{

        let deltaOD1 = log(Double(PD2_IR1) / Double(PD1_IR1));
        let deltaOD2 = log(Double(PD2_IR2) / Double(PD1_IR2));

        let changeHbRVal = (alphaHbO1 * (deltaOD2 /  Double(DPF2)) - alphaHbO2 * (deltaOD1 / Double(DPF1))) / (L2*(alphaHbR2 * alphaHbO1 - alphaHbR2 * alphaHbO2));

        return changeHbRVal;
    }
    
    public static func changeHbOHbR(_ diff_hbr_wl: Double, diff_hbo_wl: Double, distance: Double,  path: Int) -> [Double]{

        // https://mail.nmr.mgh.harvard.edu/pipermail//homer-users/2006-July/000124.html
        // https://metabolight.org/wp-content/uploads/2018/02/2014-Scholkmann-et-al.-A-review-on-continuous-wave-functional-near-infrared-spectroscopy-and-imaging-instrumentation-and-methodology.pdf
        //
        //        dOD_w1 = -log(I_w1/mean(I_w1));
        //        dOD_w2 = -log(I_w2/mean(I_w2));
        //
        //        E = GetExtinctions([850 740]);
        //        E=E(1:2,1:2);
        //
        //        dOD_w1_L = dOD_w1 * ppf/(dpf*L);
        //        dOD_w2_L = dOD_w2 * ppf/(dpf*L);
        //
        //        dOD_L = [dOD_w1_L dOD_w2_L];
        //
        //        Einv = inv(E'*E)*E';
        //
        //        HbO_HbR = Einv * dOD_L';
        //
        //        HbO = HbO_HbR(1,:);
        //        HbR = HbO_HbR(2,:);
        //
        //PD - Photodiode wavelength (e.g 740, 850, 940)
        //IR - IR Distance - (e.g. 10, 27)

        let L1 = distance; // mm

        //A = L * absorptionCoeff
        //extinction coefficients
        //https://www.researchgate.net/publication/7910860_Extinction_coefficients_of_hemoglobin_for_near-infrared_spectroscopy_of_tissue

        //Wavelength 1 = 740nm
        let alphaHbO_WL740 = 0.45;

        //Wavelength 1 = 740nm
        let alphaHbR_WL740 = 1.5;

        //Wavelength 2 = 940nm
//        let alphaHbO_WL940 = 1.35;

        //Wavelength 2 = 940nm
//        let alphaHbR_WL940 = 0.9;

        //Wavelength 3 = 850nm
        let alphaHbR_WL850 = 1.1;

        //Wavelength 3 = 850nm
        let alphaHbO_WL850 = 0.75;

        //alt 740nm        //alt 850nm         //alt 660nm           //alt 880nm            //alt 940nm
        //740 HbR = 1.5    //850 HbR = 1.1f    //660nm HbR = 1.7f    //880nm HbR = 0.8f     //940nm HbR = 0.9f
        //740 HbO = 0.45   //850 HbO = 0.75f   //660nm HbO = 0.35f   //880nm HbO = 1.25f    //940nm HbO = 1.35f

        //function for determining DPF based on gender and age
//        let age = 31;
        //DPF value
        //https://www.nature.com/articles/pr19962544/tables/3
        let DPF1 = 6.65;     //740nm
//        let DPF2 = 6.208;    //940nm
        let DPF3 = 6.12;     //850nm

        let ppf = 0.1;
        var meanhboPD1 = 0.0;
        var meanhbrPD1 = 0.0;
        
        //based on ratio to get to common 2mW/sr for LED intensity at 20mA
        let factor740nm = 3.33;
//        let factor940nm = 1.67;
        let factor850nm = 1.0;
        
        //based on photodiode used
        //https://www.mouser.ca/datasheet/2/427/vemd8081-1923765.pdf
        let qEfficiency740nm = 0.94;
        let qEfficiency850nm = 1.0;
        
        //sensitivity area for photodiode, 2x PD on long path
        let shortPathSA = 1.0;
        let longPathSA = 2.0;

        if (path == 0) {
            sumTotalHbOShort = sumTotalHbOShort + (diff_hbo_wl * (factor850nm/qEfficiency850nm))/shortPathSA;
            sumTotalHbRShort = sumTotalHbRShort + (diff_hbr_wl * (factor740nm/qEfficiency740nm))/shortPathSA;
            countTotalShort = countTotalShort + 1.0;
            meanhboPD1 = sumTotalHbOShort / countTotalShort; // hbo, > 800nm
            meanhbrPD1 = sumTotalHbRShort / countTotalShort; // hbr, < 800nm
        } else if (path == 1){
            sumTotalHbOLong = sumTotalHbOLong + (diff_hbo_wl * (factor850nm/qEfficiency850nm))/longPathSA;
            sumTotalHbRLong = sumTotalHbRLong + (diff_hbr_wl * (factor740nm/qEfficiency740nm))/longPathSA;
            countTotalLong = countTotalLong + 1.0;
            meanhboPD1 = sumTotalHbOLong / countTotalLong; // hbo, > 800nm
            meanhbrPD1 = sumTotalHbRLong / countTotalLong; // hbr, < 800nm
        }

        let E: [Double] = [
            [alphaHbO_WL850, alphaHbR_WL850],
            [alphaHbO_WL740, alphaHbR_WL740]
        ].flatMap { $0 }
        
        //creates log(n) curve with iterative mean
        var dOD_hbo_WL = 0.0
        var dOD_hbr_WL = 0.0
        
        if (path == 0) {
            dOD_hbo_WL = -1 * log(abs(diff_hbo_wl*(factor850nm/qEfficiency850nm)/shortPathSA)/meanhboPD1);
            dOD_hbr_WL = -1 * log(abs(diff_hbr_wl*(factor740nm/qEfficiency740nm)/shortPathSA)/meanhbrPD1);
        } else if (path == 1){
            dOD_hbo_WL = -1 * log(abs(diff_hbo_wl*(factor850nm/qEfficiency850nm)/longPathSA)/meanhboPD1);
            dOD_hbr_WL = -1 * log(abs(diff_hbr_wl*(factor740nm/qEfficiency740nm)/longPathSA)/meanhbrPD1);
        }

        let dOD_hbo_L1 = dOD_hbo_WL * ppf / (DPF1 * L1);
        let dOD_hbr_L1 = dOD_hbr_WL * ppf / (DPF3 * L1);

        let dOD_L: [Double] = [
            [dOD_hbo_L1, dOD_hbr_L1],
            [0.0, 0.0]
        ].flatMap { $0 }

        //Einv = inv(E'*E)*E';
        let Et = transpose2x2(E);
        let EtxE = vDSP.multiply(E, Et)
        let invEtxE = invert2x2(EtxE);
        let Einv = vDSP.multiply(invEtxE,E)
//        print(Einv)

        //HbO_HbR = Einv * dOD_L']
        var hboHbr = [0.0,0.0]

        let t_dOD_L = transpose2x1(dOD_L);
//        print(t_dOD_L)
        let HbOHbR = multiply2x2by2x1(Einv, b: t_dOD_L);
//        print(HbOHbR)
        let output = HbOHbR
        let hbo = output[0];
        let hbr = output[1];

        hboHbr[0] = hbo;
        hboHbr[1] = hbr;

        return hboHbr;
    }
    
    public static func transpose2x2(_ a: [Double]) -> [Double]{
        
        let b: [Double] = [
            [a[0],a[2]],
            [a[1],a[3]]
        ].flatMap { $0 }

        return b;
    }
    
    public static func transpose2x1(_ a: [Double]) -> [Double]{
        let b: [Double] = [
            [a[0]],
            [a[1]]
        ].flatMap { $0 }

        return b;
    }
    
    public static func invert2x2(_ a: [Double]) -> [Double]{
        var b: [Double] = [
            [0.0,0.0],
            [0.0,0.0]
        ].flatMap { $0 }

        let det1_A = 1/(a[0]*a[3] - a[1]*a[2]);

        b[0] = a[3]*det1_A
        b[1] = -1*a[1]*det1_A
        b[2] = -1*a[2]*det1_A
        b[3] = a[0]*det1_A

        return b;
    }
    
    public static func multiply2x2by2x1(_ a: [Double], b: [Double]) -> [Double]{
        var c = [0.0,0.0];

        c[0] = a[0]*b[0]+a[1]*b[1];
        c[1] = a[2]*b[0]+a[3]*b[1];

        return c;
    }
    
    public static func shortChannelCorrectionSingle(_ data10mm: Double, data27mm: Double) -> Double {

        //https://github.com/mne-tools/mne-nirs/blob/a0e7cd25f1e829826b3a9f3c3f4823dd947b5dae/mne_nirs/signal_enhancement/_short_channel_correction.py
        //        # Eqn 27 Scholkmann et al 2014
        //        alfa = np.dot(A_s, A_l) / np.dot(A_s, A_s)
        //
        //        # Eqn 26 Scholkmann et al 2014
        //        correctedVal = A_l - alfa * A_s
        //
        var alfa = 0.0;
        var A_l = 0.0;
        var A_s = 0.0;
        var correctionVal = 0.0;
        let factor = 1.000000001; //scales adjustment

        A_l = data27mm;
        A_s = data10mm;
        alfa = ((A_l * A_s) / (A_s * A_s))*factor;
        correctionVal = (A_l*factor) - (alfa * (A_s*factor));   //scale single values by factor
        correctionVal = correctionVal/factor;                   //reverse factor scaling

        return correctionVal;
    }
    
    public static func calculateDPF(_ age: Double) -> [Double]{
        var dpfVals = [0.0,0.0,0.0];

        let dpf740nm = 5.3 + pow(0.049*age, 0.723);
        let dpf940nm = 4.2 + pow((0.07*age),0.9);
        let dpf850nm = 4.5 + pow((0.058*age),0.823);

        dpfVals[0] = dpf740nm;
        dpfVals[1] = dpf940nm;
        dpfVals[2] = dpf850nm;

        return dpfVals;
    }

}
