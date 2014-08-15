package com.lomolith.common.util;

public class Converter {
    public static boolean ERROR=false;
    
    public Converter() {
    }
    
    public static double toDouble(String s) throws Exception {
        try {
            double d = Double.valueOf(s).doubleValue();
            if (Double.isNaN(d) || Double.isInfinite(d)) {
                if (ERROR) {
                    throw new ArithmeticException("Illegal double value: " + s);
                }
                else return 0.0;
            }
            return(d);
        }
        catch (Exception e) {
            if (ERROR) throw new Exception(e);
            else return 0.0;
        }
    }

    public static int toInt(Object s) throws Exception {
        try {
            return Integer.parseInt(s.toString());
        }
        catch (Exception e) {
            if (ERROR)throw new Exception(e);
            else return 0;
        }
    }
}
