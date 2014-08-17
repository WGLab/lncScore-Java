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
    Predictor p = new Predictor();
    
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

//    @Test
    public void testParseCommandLine() {
        String argv[]={"features.txt"};
        try {
            p.parse_command_line(argv);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
//    @Test
    public void testReadProblem() {
        String argv[]={"c:/workspace/ncRNA/features.txt"};
        try {
            p.parse_command_line(argv);
            p.read_problem();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    @Test
    public void testRun() {
        String argv[]={"c:/workspace/ncRNA/featurelet.txt"};
        try {
            p.run(argv);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}