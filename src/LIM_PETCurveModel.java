package pmclass.lib.pmod.models;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import pmclass.lib.pmdef.PMsettings;
import pmclass.lib.pmdef.PMunits;
import pmclass.lib.pmkin.PKgenericModel;
import pmclass.lib.pmkin.PKparameter;
import pmclass.lib.pmod.*;
import pmclass.lib.pmdata.PMdynamicData;

import cern.jet.stat.Descriptive; 
import cern.colt.list.DoubleArrayList; 

import ij.IJ;
import ij.ImagePlus;
import ij.process.*;
import ij.process.ImageProcessor.*;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import java.util.Locale;
import javax.swing.JFileChooser;

/**
 *
 * @author José María Mateos (<a href="mailto:jmmateos@mce.hggm.es">jmmateos@mce.hggm.es</a>)
 *
 * Dynamic PET segmentation algorithm based on a priori known TACs.
 * Contains many unused functions, but this is the code that was used
 * at the time of publication. No changes have been made to it, other than
 * removing commented lines.
 * 
 */
public class LIM_PETCurveModel extends PMmodel implements PMsettings {

    private final static String modelName = "LIM PET Input Function";
    private final static int PREP_N_PARAMETERS = 3;
    private final static int PREP_CHOOSE_MODEL = 0;
    private final static int PREP_CHOOSE_PATH = 1;
    private final static int PREP_CHOOSE_ITER = 2;
    private final static int PREP_N_DER_PARAMETERS = 0;

    /*
     * Parameters
     */
    private final static int N_PARAMETERS = 5;
    private final static int LV = 0;
    private final static int RV = 1;
    private final static int MC = 2;
    private final static int LVAMP = 3;
    private final static int RVAMP = 4;
    private final static int N_DER_PARAMETERS = 0;
    PMmodelData md;
    PMdynamicData dd;
    int dim_x, dim_y, dim_z, dim_t;
    int tracer_model = 0; 
    int condicion_maximo; 
    int ancho_filtro_mediana; 
    double corr_limit_lv, corr_limit_rv, corr_limit_mc; 
    double corr_factor = 1.02; 

    double lv_model[];
    double rv_model[];
    double mc_model[];
    double amp_model[];
    double lv_amp_model[];
    double rv_amp_model[];
    double mc_amp_model[];

    double iter_lv[][][];
    double iter_rv[][][];
    double iter_mc[][][];
    double amp_lv[][][];
    double amp_rv[][][];
    double amp_mc[][][];
    double iter_lv_amp[][][];
    double iter_rv_amp[][][];
    double iter_mc_amp[][][];

    public LIM_PETCurveModel() {
        super(N_PARAMETERS, N_DER_PARAMETERS, modelName);
        createPreprocesModel(PREP_N_PARAMETERS, PREP_N_DER_PARAMETERS);
        setPreprocesDefaultValues(getPreprocesModel());
    }

    /**
     *
     * @param modelData
     * @return
     * @throws PMmodelPrepareException
     */
    @Override
    public String prepare(PMmodelData modelData) throws PMmodelPrepareException {
        double [] corr_vector = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};

        md = modelData; 
        dd = modelData.dynamic; 

        long start = System.currentTimeMillis();

        tracer_model = (int) this.modelPreProc.getPar(PREP_CHOOSE_MODEL).getVal();
        int choose_file = (int) this.modelPreProc.getPar(PREP_CHOOSE_PATH).getVal();
        int max_iterations = (int) this.modelPreProc.getPar(PREP_CHOOSE_ITER).getVal();
        String path;

        if (choose_file == 1) {
            JFileChooser chooser = new JFileChooser();
            chooser.setCurrentDirectory(new File("."));
            chooser.setDialogTitle("Select directory");
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            chooser.showOpenDialog(this.getModelPreprocessingPanel());

            path = chooser.getSelectedFile().getAbsolutePath();
            path += path + System.getProperty("file.separator");

        } else {
            path = "C:/temp/";
        }

        System.out.println(path + " will be used to save the result.");

        String plasma_file = "plasma.txt";
        String tac_file = "tac.txt";

        if (tracer_model == 1) { 

            corr_limit_mc = corr_limit_lv = corr_limit_rv = 0.5;
            

        } else if (tracer_model == 2) { 

        } else if (tracer_model == 3) { 

            corr_limit_mc = 0.2;
            corr_limit_lv = 0.5;
            corr_limit_rv = 0.4;

        } else if (tracer_model == 4) { 

            corr_limit_lv = 0.5;
            corr_limit_rv = 0.1;
            corr_limit_mc = 0.5;

        } else if (tracer_model == 5) { 
            corr_limit_mc = 0.5;
            corr_limit_lv = 0.3;
            corr_limit_rv = 0;            

        } else if (tracer_model == 6) { 

            corr_limit_mc = 0.0; 
            corr_limit_lv = 0.5;
            corr_limit_rv = 0.4;
        }
        
        condicion_maximo = 0;
        ancho_filtro_mediana = 3;

        dim_x = dd.nx();
        dim_y = dd.ny();
        dim_z = dd.nz();
        dim_t = dd.ntime();

        iter_lv = new double[dim_x][dim_y][dim_z];
        iter_rv = new double[dim_x][dim_y][dim_z];
        iter_mc = new double[dim_x][dim_y][dim_z];
        amp_lv = new double[dim_x][dim_y][dim_z];
        amp_rv = new double[dim_x][dim_y][dim_z];
        amp_mc = new double[dim_x][dim_y][dim_z];
        iter_rv_amp = new double[dim_x][dim_y][dim_z];
        iter_lv_amp = new double[dim_x][dim_y][dim_z];
        iter_mc_amp = new double[dim_x][dim_y][dim_z];

        amp_model = new double[dim_t];
        lv_amp_model = new double[dim_t];
        rv_amp_model = new double[dim_t];
        mc_amp_model = new double[dim_t];


        double new_lv_model[];
        double new_rv_model[];
        double new_mc_model[];

        double temp[] = md.tBegin;


        double[][][] mask_lv;
        double[][][] mask_rv;
        double[][][] mask_mc;
        byte[][][] interior_mc = null; 

        DoubleArrayList array_lv;
        DoubleArrayList array_rv;
        DoubleArrayList array_mc;

        int it_number = 1;
        int n_lv = -1, n_rv = -1, n_mc = -1;
        int old_n_lv = -1, old_n_rv = -1, old_n_mc = -1;
        int changes = -1;
        int total = 0;

        boolean force_roi = false;
        boolean roi_has_been_forced = false;


        while ((it_number <= max_iterations)
                && (changes != 0)) {
            System.out.println("\t- Iteration " + it_number);
            this.progress.set("Iteration " + it_number);
            this.progress.setNum(100);
            this.progress.set(0);
            this.progress.setAlwaysOnTop(false);
            
            if (!roi_has_been_forced) {
                corr_limit_mc = corr_limit_lv = corr_limit_rv = (it_number > corr_vector.length)
                        ? corr_vector[corr_vector.length - 1]
                        : corr_vector[it_number - 1];
            } 


            if (it_number == 1) 
            {
                lv_model = LVCurve(temp);
                rv_model = RVCurve(temp);
                mc_model = miocardyumCurve(temp);
            } else {
                mask_lv = makeMask(iter_lv, dim_x, dim_y, dim_z);
                this.progress.set(10);
                mask_rv = makeMask(iter_rv, dim_x, dim_y, dim_z);
                this.progress.set(20);
                mask_mc = makeMask(iter_mc, dim_x, dim_y, dim_z);
                this.progress.set(30);

                if ((old_n_lv == 0 || old_n_rv == 0) && it_number == 2) {
                    force_roi = roi_has_been_forced = true;
                }

                if (force_roi) {
                    
                    System.out.println("ROI FORCED!!! SHOULD NOT HAPPEN!!!");
                    
                    corr_limit_mc = corr_limit_lv = corr_limit_rv = 0.2; // PMB
                    
                    interior_mc = getROIInteriorX(mask_mc);

                    for (int i = 0; i < dim_x; i++) {
                        for (int j = 0; j < dim_y; j++) {
                            for (int k = 0; k < dim_z; k++) {
                                if (interior_mc[i][j][k] == 0) {
                                    mask_lv[i][j][k] = 0;
                                    if (old_n_rv == 0) {

                                        mask_rv[i][j][k] = 1;
                                    }

                                } else {
                                    mask_rv[i][j][k] = 0;
                                    if (old_n_lv == 0) {

                                        mask_lv[i][j][k] = 1;
                                    }
                                }
                            }
                        }
                    }
                    force_roi = false;
                }

                this.progress.set(40);


                new_lv_model = new double[dim_t];
                new_rv_model = new double[dim_t];
                new_mc_model = new double[dim_t];

                for (int i = 0; i < dim_x; i++) {
                    for (int j = 0; j < dim_y; j++) {
                        for (int k = 0; k < dim_z; k++) {
                            if (mask_rv[i][j][k] > 0) {
                                n_rv++;                                
                                for (int t = 0; t < dim_t; t++) {
                                    new_rv_model[t] += dd.getV(i, j, k, t);
                                }
                            }

                            if (mask_lv[i][j][k] > 0) {
                                n_lv++;
                                for (int t = 0; t < dim_t; t++) {
                                    new_lv_model[t] += dd.getV(i, j, k, t);
                                }
                            }

                            if (mask_mc[i][j][k] > 0) {
                                n_mc++;
                                for (int t = 0; t < dim_t; t++) {
                                    new_mc_model[t] += dd.getV(i, j, k, t);
                                }
                            }
                        }
                    }
                }

                if (roi_has_been_forced) {


                    int max_rv = getPositionMax(new_rv_model);
                    int max_lv = getPositionMax(new_lv_model);

                    if (max_rv > max_lv) {
                        double aux[] = new_rv_model;
                        new_rv_model = new_lv_model;
                        new_lv_model = aux;
                    }
                }

                System.out.println("\t- Pixels (LV, RV, MC): "
                        + n_lv + ", " + n_rv + ", " + n_mc);

                for (int t = 0; t < dim_t; t++) {
                    if (n_rv > 0) {
                        rv_model[t] = new_rv_model[t] / n_rv;
                    }
                    if (n_lv > 0) {
                        lv_model[t] = new_lv_model[t] / n_lv;
                    }
                    if (n_mc > 0) {
                        mc_model[t] = new_mc_model[t] / n_mc;
                    }
                }

            }


            this.progress.set(50);

            if (it_number != 1) {
                changes = Math.abs(old_n_rv - n_rv) + Math.abs(old_n_lv - n_lv)
                        + Math.abs(old_n_mc - n_mc);
                if ((old_n_rv == 0 || old_n_lv == 0 || old_n_mc == 0) && changes == 0) {
                    changes = 1;
                }
            }
            total = old_n_rv + old_n_lv + old_n_mc;
            
            old_n_rv = n_rv;
            old_n_lv = n_lv;
            old_n_mc = n_mc;
            n_rv = n_lv = n_mc = 0;

            array_lv = new DoubleArrayList(lv_model);
            array_rv = new DoubleArrayList(rv_model);
            array_mc = new DoubleArrayList(mc_model);

            itera(array_lv, array_rv, array_mc);
            it_number++; 

        }

        System.out.println("\t- Iteration converged");

        double max_lv = 0.0, max_rv = 0.0, max_mc = 0.0;

        for (int i = 0; i < dim_x; i++) {
            for (int j = 0; j < dim_y; j++) {
                for (int k = 0; k < dim_z; k++) {
                    if (iter_lv[i][j][k] > max_lv) {
                        max_lv = iter_lv[i][j][k];
                    }
                    if (iter_rv[i][j][k] > max_rv) {
                        max_rv = iter_rv[i][j][k];
                    }
                    if (iter_mc[i][j][k] > max_mc) {
                        max_mc = iter_mc[i][j][k];
                    }
                }
            }
        }

        double final_lv[], final_rv[], final_mc[];
        n_lv = n_rv = n_mc = 0;
        final_lv = new double[dim_t];
        final_rv = new double[dim_t];
        final_mc = new double[dim_t];

        double sel_limit_lv = max_lv * 0.7;
        double sel_limit_rv = max_rv * 0.7;
        double sel_limit_mc = max_mc * 0.7;


        for (int i = 0; i < dim_x; i++) {
            for (int j = 0; j < dim_y; j++) {
                for (int k = 0; k < dim_z; k++) {
                    if (iter_lv[i][j][k] > sel_limit_lv) {
                        n_lv++;
                        for (int t = 0; t < dim_t; t++) {
                            final_lv[t] += dd.getV(i, j, k, t);
                        }
                    } else if (iter_rv[i][j][k] > sel_limit_rv) {
                        n_rv++;
                        for (int t = 0; t < dim_t; t++) {
                            final_rv[t] += dd.getV(i, j, k, t);
                        }
                    } else if (iter_mc[i][j][k] > sel_limit_mc) {
                        n_mc++;
                        for (int t = 0; t < dim_t; t++) {
                            final_mc[t] += dd.getV(i, j, k, t);
                        }
                    }
                }
            }
        }

        for (int t = 0; t < dim_t; t++) {
            final_lv[t] /= n_lv;
            final_rv[t] /= n_rv;
            final_mc[t] /= n_mc;
        }

        for (int i = 0; i < dim_x; i++) {
            for (int j = 0; j < dim_y; j++) {
                for (int k = 0; k < dim_z; k++) {
                    for (int t = 0; t < dim_t; t++) {
                        amp_model[t] = dd.getV(i, j, k, t);
                    }
                    if (iter_lv[i][j][k] > sel_limit_lv) {
                        amp_lv[i][j][k] = getMax(amp_model);
                    }
                    if (iter_rv[i][j][k] > sel_limit_rv) {
                        amp_rv[i][j][k] = getMax(amp_model);
                    }
                    if (iter_mc[i][j][k] > sel_limit_mc) {
                        amp_mc[i][j][k] = getMax(amp_model);
                    }
                }
            }
        }

        double max_amp_lv, max_amp_rv;
        double max_amp_mc;

        max_amp_lv = getMaxVol(amp_lv, dim_x, dim_y, dim_z);
        max_amp_rv = getMaxVol(amp_rv, dim_x, dim_y, dim_z);
        max_amp_mc = getMaxVol(amp_mc, dim_x, dim_y, dim_z);

        double amp_factor = 0.5;

        int lv_amp_n = 0, rv_amp_n = 0;
        int mc_amp_n = 0;

        for (int i = 0; i < dim_x; i++) {
            for (int j = 0; j < dim_y; j++) {
                for (int k = 0; k < dim_z; k++) {
                    if (amp_lv[i][j][k] >= max_amp_lv * amp_factor) {
                        lv_amp_n++;
                        for (int t = 0; t < dim_t; t++) {
                            this.lv_amp_model[t] += dd.getV(i, j, k, t);
                        }

                    } else {
                        amp_lv[i][j][k] = 0;
                    }

                    if (amp_rv[i][j][k] >= max_amp_rv * amp_factor) {
                        rv_amp_n++;
                        for (int t = 0; t < dim_t; t++) {
                            this.rv_amp_model[t] += dd.getV(i, j, k, t);
                        }

                    } else {
                        amp_rv[i][j][k] = 0;
                    }
                    if (amp_mc[i][j][k] >= max_amp_mc * amp_factor) {
                        mc_amp_n++;
                        for (int t = 0; t < dim_t; t++) {
                            this.mc_amp_model[t] += dd.getV(i, j, k, t);
                        }

                    } else {
                        amp_mc[i][j][k] = 0;
                    }
                    
                }
            }
        }
        for (int t = 0; t < dim_t; t++) {
            lv_amp_model[t] /= lv_amp_n;
            rv_amp_model[t] /= rv_amp_n;
            mc_amp_model[t] /= rv_amp_n;
        }


        DecimalFormat df = new DecimalFormat("######", new DecimalFormatSymbols(Locale.ENGLISH));
        df.setMaximumFractionDigits(6);
        df.setMaximumIntegerDigits(6);
        df.setMinimumFractionDigits(0);
        df.setMinimumIntegerDigits(1);
        df.setParseIntegerOnly(false);
        df.setDecimalSeparatorAlwaysShown(true);

        System.out.println("Result curves are being stored on " + path + plasma_file + " and "
                + path + tac_file + "\nusing the following correlations as the lower limit:");
        System.out.println("LV: " + df.format(sel_limit_lv));
        System.out.println("RV: " + df.format(sel_limit_rv));
        System.out.println("MC: " + df.format(sel_limit_mc));


        try {
            FileWriter plasma_wr = new FileWriter(path + plasma_file);
            plasma_wr.write("time[seconds]\tvalue[Bq/ml]\n");
            for (int t = 0; t < dim_t; t++) {
                String stemp = "\t" + df.format(md.tBegin[t] + (md.tEnd[t] - md.tBegin[t]) / 2.0) + "\t"
                        + df.format(final_lv[t]) + "\n";
                plasma_wr.write(stemp);
            }
            plasma_wr.close();
        } catch (IOException ex) {
            Logger.getLogger(LIM_PETCurveModel.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            FileWriter tac_wr = new FileWriter(path + tac_file);
            tac_wr.write("start[seconds] \tend[kBq/cc]\tLV\tRV\tMC\tLVamp\tRVamp\tMCamp\n");
            for (int t = 0; t < dim_t; t++) {
                tac_wr.write(df.format(md.tBegin[t]) + "\t"
                        + df.format(md.tEnd[t]) + "\t"
                        + df.format(final_lv[t]) + "\t"
                        + df.format(final_rv[t]) + "\t"
                        + df.format(final_mc[t]) + "\t"
                        + df.format(lv_amp_model[t]) + "\t"
                        + df.format(rv_amp_model[t]) + "\t"
                        + df.format(mc_amp_model[t]) + "\n");
            }
            tac_wr.close();
        } catch (IOException ex) {
            Logger.getLogger(LIM_PETCurveModel.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println("Time elapsed: " + (double) (System.currentTimeMillis() - start) / 1000
                + " seconds");

        return "Done";

    }

    /**
     *
     * @param timeVect
     * @param paramVect
     * @param x
     * @param y
     * @param z
     * @throws PMmodelCalcException
     */
    @Override
    public void calculate(double[] timeVect, double[] paramVect, int x, int y, int z) throws PMmodelCalcException {

        paramVect[LV] = iter_lv[x][y][z];
        paramVect[RV] = iter_rv[x][y][z];
        paramVect[MC] = iter_mc[x][y][z];
        paramVect[LVAMP] = amp_lv[x][y][z];
        paramVect[RVAMP] = amp_rv[x][y][z];

    }

    /**
     *
     */
    @Override
    public void finalizeModel() {
    }

    /**
     *
     * @param modelData
     * @throws PMmodelSetDataException
     */
    @Override
    public void setModelData(PMmodelData modelData) throws PMmodelSetDataException {
    }

    /**
     *
     * @return
     */
    @Override
    public PMmodelData getModelData() {
        return this.md;
    }

    /**
     *
     * @return
     */
    @Override
    public int[] getOptionSettings() {
        int i;
        int options[] = new int[PMOD_NR_OPTIONS];
        for (i = 0; i < PMOD_NR_OPTIONS; i++) //disable all options; are selectivel
        {
            options[i] = PMOD_OPT_DISABLED;
        }

        options[PMOD_OPT_MDL_PREPROC] = PMOD_OPT_MANDATORY;

        return (options);

    }

    /**
     *
     * @return
     */
    @Override
    public String[] getOptionNames() {
        return null;
    }

    /**
     *
     * @return
     */
    @Override
    public boolean canBeParallel() {
        return false;
    }

    /**
     *
     * @param model
     */
    @Override
    protected void setDefaultValues(PKgenericModel model) {
        PKparameter par;

        par = model.getPar(LV);
        par.setName("Left Ventricle");
        par.setUnit(PMunits.VU_1_1);
        par.setTooltip("LV");
        par.setFit(true);

        par = model.getPar(RV);
        par.setName("Right Ventricle");
        par.setUnit(PMunits.VU_1_1);
        par.setTooltip("RV");
        par.setFit(true);

        par = model.getPar(MC);
        par.setName("Myocardium");
        par.setUnit(PMunits.VU_1_1);
        par.setTooltip("MC");
        par.setFit(true);

        par = model.getPar(LVAMP);
        par.setName("Left ventricle - amplitude corrected");
        par.setUnit(PMunits.VU_1_1);
        par.setTooltip("LV Amp");
        par.setFit(true);

        par = model.getPar(RVAMP);
        par.setName("Right ventricle - amplitude corrected");
        par.setUnit(PMunits.VU_1_1);
        par.setTooltip("RV Amp");
        par.setFit(true);
    }

    private void setPreprocesDefaultValues(PKgenericModel preprocesModel) {
        PKparameter par;

        par = preprocesModel.getPar(PREP_CHOOSE_MODEL);
        String instructions = "Choose model to use \n(1 - NH3, 2/3 - Rb/Noisy Rb, 4 - FDG, 5 - FDHR"
                + ", 6 - O15)";
        par.setName(instructions);
        par.setUnit(PMunits.UNDEF);
        par.setTooltip(instructions);
        par.setVal(1.0);
        par.setFit(true);
        par.setFitable(false);
        par.setRestrict(false);

        par = preprocesModel.getPar(PREP_CHOOSE_PATH);
        instructions = "Set to 1 to choose a directory for saving data or 0 to save in C:/temp/";
        par.setName(instructions);
        par.setUnit(PMunits.UNDEF);
        par.setTooltip(instructions);
        par.setVal(0.0);
        par.setFit(true);
        par.setFitable(false);
        par.setRestrict(false);

        par = preprocesModel.getPar(PREP_CHOOSE_ITER);
        instructions = "Set iterations number (default = 9)";
        par.setName(instructions);
        par.setUnit(PMunits.UNDEF);
        par.setTooltip(instructions);
        par.setVal(9.0);
        par.setFit(true);
        par.setFitable(false);
        par.setRestrict(false);


    }

    private double[] normalize(double[] data) {
        double result[] = new double[data.length];
        double max = 0.0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] > max) {
                max = data[i];
            }
        }

        for (int i = 0; i < data.length; i++) {
            result[i] = data[i] / max;
        }

        return result;
    }

    private double getMax(double[] data) {
        double max = -20.0; 
        for (int i = 0; i < data.length; i++) {
            if (data[i] > max) {
                max = data[i];
            }
        }
        return max;

    }

    private double getMaxVol(double[][][] data, int i, int j, int k) {
        double max = -20.0;
        for (int x = 0; x < i; x++) {
            for (int y = 0; y < j; y++) {
                for (int z = 0; z < k; z++) {
                    if (data[x][y][z] > max) {
                        max = data[x][y][z];
                    }
                }
            }
        }
        return max;
    }

    private int getPositionMax(double[] data) {

        double max = -20.0; 
        int index = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] > max) {
                max = data[i];
                index = i;
            }
        }
        return index;


    }

    private double getMaxIndex(double[] data, int from, int to) {
        double max = -20.0; 
        for (int i = from; i < to; i++) {
            if (data[i] > max) {
                max = data[i];
            }
        }
        return max;

    }

    private double getMin(double[] data) {
        double min = 200000;
        for (int i = 0; i < data.length; i++) {
            if (data[i] < min) {
                min = data[i];
            }
        }
        return min;

    }

    private double abs(double data) {
        if (data >= 0) {
            return data;
        } else {
            return -data;
        }
    }

    private double truncate(double data) {
        if (data >= 0) {
            return data;
        } else {
            return 0;
        }
    }

    private double[] miocardyumCurve(double temp[]) {
        double resultado[] = new double[temp.length];        
        
        
        double a = 202.59315269, b = 41.17458748, c = 0.09587109;

        for (int i = 0; i < resultado.length; i++) {
            resultado[i] = a * (1 - Math.exp(-temp[i] / b)) + c * temp[i]; 
        }

        return resultado;

    }

    private double[] LVCurve(double temp[]) {
        double resultado[] = new double[temp.length];
        
        
        double a = 0.04404489, b = 3.50217765, c = 10.14307199;
        
        for (int i = 0; i < resultado.length; i++) {
            resultado[i] = a * Math.pow(temp[i], b) * Math.exp(-temp[i] / c); 
        }

        return resultado;

    }

    private double[] RVCurve(double temp[]) {
        double resultado[] = new double[temp.length];
                
        double a = 0.9944634, b = 3.6728159, c = 4.7995968;
        
        
        
        for (int i = 0; i < resultado.length; i++) {
            resultado[i] = a * Math.pow(temp[i], b) * Math.exp(-temp[i] / c); 
        }

        return resultado;

    }

    private void itera(DoubleArrayList array_lv, DoubleArrayList array_rv,
            DoubleArrayList array_mc) {

        double dev_lv = Descriptive.standardDeviation(
                Descriptive.sampleVariance(array_lv, Descriptive.mean(array_lv)));
        double dev_rv = Descriptive.standardDeviation(
                Descriptive.sampleVariance(array_rv, Descriptive.mean(array_rv)));
        double dev_mc = Descriptive.standardDeviation(
                Descriptive.sampleVariance(array_mc, Descriptive.mean(array_mc)));

        double timeVect[] = new double[dim_t];

        for (int i = 0; i < dim_x; i++) {
            this.progress.set((int) (((double) i / (double) dim_x) * 50.0 + 50.0));
            for (int j = 0; j < dim_y; j++) {
                for (int k = 0; k < dim_z; k++) {

                    for (int t = 0; t < dim_t; t++) {
                        timeVect[t] = dd.getV(i, j, k, t);
                    }


                    if ((abs(getMin(timeVect)) < getMax(timeVect))
                            && (getMax(timeVect) > condicion_maximo)) {

                        DoubleArrayList array_tv = new DoubleArrayList(normalize(smooth(timeVect)));
                        double dev_tv = Descriptive.standardDeviation(
                                Descriptive.sampleVariance(array_tv, Descriptive.mean(array_tv)));
                        double corr_lv = Descriptive.correlation(array_lv, dev_lv, array_tv, dev_tv);
                        double corr_rv = Descriptive.correlation(array_rv, dev_rv, array_tv, dev_tv);
                        double corr_mc = Descriptive.correlation(array_mc, dev_mc, array_tv, dev_tv);

                        double values[] = new double[]{corr_lv, corr_rv, corr_mc};
                        double aux = getMax(values);

                        if (tracer_model != 2 && tracer_model != 3) {

                            if ((aux == corr_lv) && aux > corr_limit_lv) {
                                iter_lv[i][j][k] = aux;
                                iter_rv[i][j][k] = iter_mc[i][j][k] = 0;
                            } else if (aux == corr_rv && aux > corr_limit_rv) {
                                iter_rv[i][j][k] = aux;
                                iter_lv[i][j][k] = iter_mc[i][j][k] = 0;
                            } else if (aux == corr_mc && aux > corr_limit_mc) {
                                iter_mc[i][j][k] = aux;
                                iter_rv[i][j][k] = iter_lv[i][j][k] = 0;
                            } else {
                                iter_rv[i][j][k] = iter_lv[i][j][k] = iter_rv[i][j][k] = 0;
                            }

                        } else { 

                            double amp = getMax(timeVect);
                            double rel = amp / Math.abs(timeVect[timeVect.length - 1]);

                            if ((aux == corr_lv || aux == corr_mc) && aux > corr_limit_lv && rel > 6) {
                                iter_lv[i][j][k] = aux;
                                iter_rv[i][j][k] = iter_mc[i][j][k] = 0;
                            } else if (aux == corr_rv && aux > corr_limit_rv) {
                                iter_rv[i][j][k] = aux;
                                iter_lv[i][j][k] = iter_mc[i][j][k] = 0;
                            } else if ((aux == corr_mc || aux == corr_lv) && aux > corr_limit_mc && rel < 4) {
                                iter_mc[i][j][k] = aux;
                                iter_rv[i][j][k] = iter_lv[i][j][k] = 0;
                            } else {
                                iter_rv[i][j][k] = iter_lv[i][j][k] = iter_rv[i][j][k] = 0;
                            }

                        }
                    }


                }
            }
        }

    }

    private double[][][] makeMask(double[][][] data, int dim_x, int dim_y, int dim_z) {

        double[][][] resultado = new double[dim_x][dim_y][dim_z];
        float rodaja[][] = new float[dim_x][dim_y];
        ImageProcessor ip;
        ImagePlus image;
        for (int i = 0; i < dim_z; i++) {

            for (int j = 0; j < dim_x; j++) {
                for (int k = 0; k < dim_y; k++) {
                    rodaja[j][k] = (float) data[j][k][i];
                }
            }
            ip = new FloatProcessor(rodaja);
            image = new ImagePlus("data", ip);

            
            IJ.run(image, "Median...", "radius=" + ancho_filtro_mediana);
            ip.threshold(0);

            for (int j = 0; j < dim_x; j++) {
                for (int k = 0; k < dim_x; k++) {
                    resultado[j][k][i] = ip.getPixel(j, k);
                }
            }

        }

        return resultado;
    }

    private byte[][][] getROIInteriorX(double[][][] mask_mc) {

        byte[][][] area = new byte[dim_x][dim_y][dim_z];
        double temp[];
        int coincidencial = 0, coincidenciar = 0, coincidenciat = 0, coincidenciab = 0;

        for (int rodaja = 1; rodaja < dim_x - 1; rodaja++) // me quito los bordes
        {
            for (int i = 1; i < dim_z - 1; i++) {
                for (int j = 1; j < dim_y - 1; j++) {
                    if ((int) mask_mc[rodaja][j][i] == 0) { 

                        coincidencial = coincidenciar = coincidenciat = coincidenciab = 0;

                        temp = getRow(mask_mc, dim_z, j, rodaja);
                        if (getMaxIndex(temp, 0, i) > 0) {
                            coincidencial++;
                        }
                        if (getMaxIndex(temp, i + 1, temp.length) > 0) {
                            coincidenciar++;
                        }

                        temp = getColumn(mask_mc, i, dim_y, rodaja);
                        if (getMaxIndex(temp, 0, j) > 0) {
                            coincidenciat++;
                        }
                        if (getMaxIndex(temp, j + 1, temp.length) > 0) {
                            coincidenciab++;
                        }

                        area[rodaja][j][i] = (coincidencial + coincidenciar
                                + coincidenciat + coincidenciab >= 3)
                                ? (byte) 1 : (byte) 0;

                    }

                }
            }
        }
        return area;
    }

    private double[] getColumn(double[][][] data, int x, int size_y, int rodaja) {

        double result[] = new double[size_y];
        for (int i = 0; i < size_y; i++) {
            result[i] = data[rodaja][i][x];
        }

        return result;

    }

    private double[] getRow(double[][][] data, int size_x, int y, int rodaja) {

        double result[] = new double[size_x];
        for (int i = 0; i < size_x; i++) {
            result[i] = data[rodaja][y][i];
        }

        return result;

    }
    
    /**
     * Smooths the data before doing any further operation with it
     * @param original Original data
     * @return Smoothed data
     */
    private double[] smooth(double[] original) {
        return smooth3(original);
    }
    
    private double[] smooth5(double[] original) {

        double result[] = original.clone();
        double kernel[] = {0.13, 0.185, 0.37, 0.185, 0.13};
        int t = result.length;
        for (int i = 2; i < t - 2; i++) {
            result[i] = result[i - 2] * kernel[0] + result[i - 1] * kernel[1] + result[i] * kernel[2]
                    + result[i + 1] * kernel[3] + result[i + 2] * kernel[4];
        }
        return result;
    }
    
    private double[] smooth3(double[] original) {

        double result[] = original.clone();
        double kernel[] = {0.25, 0.5, 0.25};
        int t = result.length;
        for (int i = 1; i < t - 1; i++) {
            result[i] = result[i - 1] * kernel[0] + result[i] * kernel[1] + 
                    result[i + 1] * kernel[2];                    
        }
        return result;
    }
}
