package com.lomolith.common.util;

public class EntropyUtil {
        private double variance;
        
        public EntropyUtil() {
        }
        
        public double getSTDDev(double[] observations) {
            double sumSqs = 0;
            double mean = 0;
            for (int i=0; i<observations.length; i++) mean += observations[i];
            mean /= (double) observations.length;
            for (int i=0; i<observations.length; i++) sumSqs += (observations[i] - mean) * (observations[i] - mean);
            variance = sumSqs / (double) (observations.length - 1);
            return Math.sqrt(variance);
	}    

        /**
         * <p>The entropy for a Gaussian-distribution random variable with
         *  variance \sigma is 0.5*\log_e{2*pi*e*\sigma}.</p>
         * 
         * <p>Here we compute the entropy assuming that the recorded estimation of the
         *  variance is correct (i.e. we will not make a bias correction for limited
         *  observations here).</p>
         * 
         * @return the entropy of the previously provided observations or from the supplied
         *   covariance matrix. Entropy returned in nats, not bits!
         */
        public double getNormDistDiffEntropy(double[] observations) {
                return 0.5 * Math.log(2.0*Math.PI*Math.E*getSTDDev(observations));
        }
}
