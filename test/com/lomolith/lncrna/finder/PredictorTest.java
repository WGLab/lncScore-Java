package com.lomolith.lncrna.finder;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author mittens
 */
public class PredictorTest {
    
    public PredictorTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void test() {
        Predictor p = new Predictor();
        String argv[]={"-pk"};
        try {
            p.parse_command_line(argv);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}