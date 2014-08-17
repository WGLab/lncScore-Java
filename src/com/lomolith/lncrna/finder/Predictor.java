package com.lomolith.lncrna.finder;

import com.lomolith.common.util.Converter;
import java.io.*;
import java.util.*;
import libsvm.*;

public class Predictor {
    private svm_parameter param;		// set by parse_command_line
    private svm_problem prob;                   // set by read_problem
    private svm_model model;
    private String input_file_name;		// set by parse_command_line
    private String model_file_name;		// set by parse_command_line
    private String error_msg;
    private int cross_validation;
    private int nr_fold;
    
    public static svm_print_interface svm_print_null = new svm_print_interface() {
        public void print(String s) {}
    };

    public static void exit_with_help() {
        System.out.print(
         "Usage: svm_train [options] training_set_file [model_file]\n"
        +"options:\n"
        +"-s svm_type : set type of SVM (default 0)\n"
        +"	0 -- C-SVC		(multi-class classification)\n"
        +"	1 -- nu-SVC		(multi-class classification)\n"
        +"	2 -- one-class SVM\n"
        +"	3 -- epsilon-SVR	(regression)\n"
        +"	4 -- nu-SVR		(regression)\n"
        +"-t kernel_type : set type of kernel function (default 2)\n"
        +"	0 -- linear: u'*v\n"
        +"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
        +"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
        +"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
        +"	4 -- precomputed kernel (kernel values in training_set_file)\n"
        +"-d degree : set degree in kernel function (default 3)\n"
        +"-g gamma : set gamma in kernel function (default 1/num_features)\n"
        +"-r coef0 : set coef0 in kernel function (default 0)\n"
        +"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
        +"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
        +"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
        +"-m cachesize : set cache memory size in MB (default 100)\n"
        +"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
        +"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
        +"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
        +"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
        +"-v n : n-fold cross validation mode\n"
        +"-q : quiet mode (no outputs)\n"
        );
        System.exit(1);
    }

    public void do_cross_validation() {
        int i;
        int total_correct = 0;
        double total_error = 0;
        double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
        double[] target = new double[prob.l];

        svm.svm_cross_validation(prob,param,nr_fold,target);
        if( param.svm_type == svm_parameter.EPSILON_SVR ||
            param.svm_type == svm_parameter.NU_SVR){
            for(i=0;i<prob.l;i++) {
                double y = prob.y[i];
                double v = target[i];
                total_error += (v-y)*(v-y);
                sumv += v;
                sumy += y;
                sumvv += v*v;
                sumyy += y*y;
                sumvy += v*y;
            }
            System.out.print("Cross Validation Mean squared error = "+total_error/prob.l+"\n");
            System.out.print("Cross Validation Squared correlation coefficient = "+
                ((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
                ((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))+"\n" );
        }
        else {
            for(i=0;i<prob.l;i++)
                if(target[i] == prob.y[i]) ++total_correct;
                System.out.print("Cross Validation Accuracy = "+100.0*total_correct/prob.l+"%\n");
        }
    }
	
    public void run(String argv[]) throws Exception {
        parse_command_line(argv);
        read_problem();
        error_msg = svm.svm_check_parameter(prob,param);
        if(error_msg != null){
            System.err.print("ERROR: "+error_msg+"\n");
            System.exit(1);
        }

        if(cross_validation != 0) do_cross_validation();
        else {
            model = svm.svm_train(prob,param);
System.out.println(model)            ;
            svm.svm_save_model(model_file_name,model);
        }
    }

    public void parse_command_line(String argv[]) throws Exception {
        int i;
        svm_print_interface print_func = null;	// default printing to stdout

        param = new svm_parameter();
        // default values
        param.svm_type = svm_parameter.C_SVC;
        param.kernel_type = svm_parameter.RBF;
        param.degree = 3;
        param.gamma = 0;	// 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 1;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = new int[0];
        param.weight = new double[0];
        cross_validation = 0;
        
        // parse options
        for(i=0;i<argv.length;i++) {
            if(argv[i].charAt(0) != '-') break;
            if(++i>=argv.length) exit_with_help();
            switch(argv[i-1].charAt(1)) {
                case 's':
                    param.svm_type = Converter.toInt(argv[i]);
                    break;
                case 't':
                    param.kernel_type = Converter.toInt(argv[i]);
                    break;
                case 'd':
                    param.degree = Converter.toInt(argv[i]);
                    break;
                case 'g':
                    param.gamma = Converter.toDouble(argv[i]);
                    break;
                case 'r':
                    param.coef0 = Converter.toDouble(argv[i]);
                    break;
                case 'n':
                    param.nu = Converter.toDouble(argv[i]);
                    break;
                case 'm':
                    param.cache_size = Converter.toDouble(argv[i]);
                    break;
                case 'c':
                    param.C = Converter.toDouble(argv[i]);
                    break;
                case 'e':
                    param.eps = Converter.toDouble(argv[i]);
                    break;
                case 'p':
                    param.p = Converter.toDouble(argv[i]);
                    break;
                case 'h':
                    param.shrinking = Converter.toInt(argv[i]);
                    break;
                case 'b':
                    param.probability = Converter.toInt(argv[i]);
                    break;
                case 'q':
                    print_func = svm_print_null;
                    i--;
                    break;
                case 'v':
                    cross_validation = 1;
                    nr_fold = Converter.toInt(argv[i]);
                    if(nr_fold < 2) {
                        System.err.print("n-fold cross validation: n must >= 2\n");
                        exit_with_help();
                    }
                    break;
                case 'w':
                    ++param.nr_weight;
                    {
                        int[] old = param.weight_label;
                        param.weight_label = new int[param.nr_weight];
                        System.arraycopy(old,0,param.weight_label,0,param.nr_weight-1);
                    }

                    {
                        double[] old = param.weight;
                        param.weight = new double[param.nr_weight];
                        System.arraycopy(old,0,param.weight,0,param.nr_weight-1);
                    }

                    param.weight_label[param.nr_weight-1] = Converter.toInt(argv[i-1].substring(2));
                    param.weight[param.nr_weight-1] = Converter.toDouble(argv[i]);
                    break;
                default:
                    System.err.print("Unknown option: " + argv[i-1] + "\n");
                    exit_with_help();
            }
        }

        svm.svm_set_print_string_function(print_func);

        // determine filenames

        if(i>=argv.length) exit_with_help();

        input_file_name = argv[i];

        if(i<argv.length-1) model_file_name = argv[i+1];
        else {
            int p = argv[i].lastIndexOf('/');
            ++p;	// whew...
            model_file_name = argv[i].substring(p)+".model";
        }
    }

    // read in a problem (in svmlight format)

    public void read_problem() throws Exception {        
        BufferedReader fp = new BufferedReader(new FileReader(input_file_name));
        Vector<Double> vy = new Vector<Double>();
        Vector<svm_node[]> vx = new Vector<svm_node[]>();

        String line=fp.readLine();
        String idx[]=line.split("\t");
        int max_index = idx.length;
        
        while((line=fp.readLine())!=null) {
            String st[]=line.split("\t");
//            StringTokenizer st = new StringTokenizer(line," \t\n\r\f:");
//            vy.addElement(Converter.toDouble(st.nextToken()));
            vy.addElement(Converter.toDouble(st[2]));
            svm_node[] x = new svm_node[idx.length+1];
            for(int j=3;j<idx.length;j++) {
                x[j] = new svm_node();
                x[j].index = Converter.toInt(j);
                x[j].value = Converter.toDouble(st[j]);
            }
//            x[idx.length] = new svm_node();
//            x[idx.length].index=-1;
            vx.addElement(x);
        }

        prob = new svm_problem();
        prob.l = vy.size();
        prob.x = new svm_node[prob.l][];
        for(int i=0;i<prob.l;i++) prob.x[i] = vx.elementAt(i);
        prob.y = new double[prob.l];
        for(int i=0;i<prob.l;i++) prob.y[i] = vy.elementAt(i);

        if(param.gamma == 0 && max_index > 0) param.gamma = 1.0/max_index;

        if(param.kernel_type == svm_parameter.PRECOMPUTED)
            for(int i=0;i<prob.l;i++) {
                if (prob.x[i][0].index != 0) {
                    System.err.print("Wrong kernel matrix: first column must be 0:sample_serial_number\n");
                    System.exit(1);
                }
                if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index) {
                    System.err.print("Wrong input format: sample_serial_number out of range\n");
                    System.exit(1);
                }
            }
            fp.close();
    }
}