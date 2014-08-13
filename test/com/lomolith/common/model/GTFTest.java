/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.common.model;

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
public class GTFTest {
    
    GTF g;
    public GTFTest() {
        g = new GTF("");
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

    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}
    @Test
    public void testIsExist() {
        
    }
    @Test
    public void testIsGTF() {
        
    }
    @Test
    public void testInitiateGTF() {
        
    }
    @Test 
    public void testFileName() {
        g.setFile("Test");
        assertEquals(g.getFile(),"Test");
    }
}
