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
    public void testExit_with_help() {
        System.out.println("exit_with_help");
        Predictor.exit_with_help();
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testDo_cross_validation() {
        System.out.println("do_cross_validation");
        Predictor instance = new Predictor();
        instance.do_cross_validation();
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testRun() throws Exception {
        System.out.println("run");
        String[] argv = null;
        Predictor instance = new Predictor();
        instance.run(argv);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testMain() throws Exception {
        System.out.println("main");
        String[] argv = null;
        Predictor.main(argv);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testAtof() {
        System.out.println("atof");
        String s = "";
        double expResult = 0.0;
        double result = Predictor.atof(s);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testAtoi() {
        System.out.println("atoi");
        String s = "";
        int expResult = 0;
        int result = Predictor.atoi(s);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testParse_command_line() {
        System.out.println("parse_command_line");
        String[] argv = null;
        Predictor instance = new Predictor();
        instance.parse_command_line(argv);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    @Test
    public void testRead_problem() throws Exception {
        System.out.println("read_problem");
        Predictor instance = new Predictor();
        instance.read_problem();
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
