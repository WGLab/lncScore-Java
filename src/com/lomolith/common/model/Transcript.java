/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.common.model;

/**
 *
 * @author mittens
 */
public class Transcript {
    public String id;
    public String gene;
    public float FPKM;
    public float conf_lo;
    public float conf_hi;
    public String method;
    public String type;
    public String strand;
    public String chr;
    public long start;
    public long end;
    public int exon;
    public String sequence;
    public String orf_seq;
    public String orf_open=null; // S: Start, E: End, P: Paired null: No detect
}
