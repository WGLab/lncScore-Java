/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.common.util;

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
public class DNAUtilTest {
    
    public DNAUtilTest() {
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
    public void testGetSeq() {
        DNAUtil du = new DNAUtil();
        try {
            String seq=du.getSeq("chr1:100000-110000");
            System.out.println(seq.length()+"/"+seq+"/");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
