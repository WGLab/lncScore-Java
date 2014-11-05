/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.annotation;

import com.lomolith.common.model.GTF;
import com.lomolith.common.model.Transcript;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mittens
 */
public class GTFManager {
    private String filename;
    
    public boolean VERBOSE=false;
    public boolean FILTER_SIZE=true;
    public int TOTAL=0;
    
    public void setFile(String file) {
        filename = file;
        
    }
    
    public GTF read() throws FileNotFoundException, IOException {
        return read(null);
    }
    
    public GTF read(GTF annot) throws FileNotFoundException, IOException {
        GTF gtf = new GTF(filename);
        gtf.FILTER_SIZE=FILTER_SIZE;
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        String line="";
        int cnt=0;        
boolean diff=false;        
        while ((line=br.readLine())!=null) {
            if (!line.startsWith("#")) {
                gtf.setTranscript(line);
                cnt++;
if (cnt!=gtf.transcripts.size()&&!diff) {System.out.println(cnt+"/"+gtf.transcripts.size()); diff=true;}
            }
        }
        gtf.sortTranscripts();
        TOTAL=cnt;
        br.close();
        fr.close();
        
        if (annot!=null) {
            System.out.print("Removing known transcripts from "+gtf.getSize()+"...");
            gtf=exclude(gtf, annot);
            System.out.println("done. "+gtf.getSize()+" transcripts remain");
        }
        return gtf;
    }
    
    private int getChr(String chr) {
        try { 
            if (chr.indexOf("_")!=-1) return Integer.parseInt(chr.substring(0,chr.indexOf("_")));
            else return Integer.parseInt(chr); 
        } catch(NumberFormatException e) { 
            if (chr.toUpperCase().startsWith("X")) return 23;
            else if (chr.toUpperCase().startsWith("Y")) return 24;
            else return 999; 
        }        
    }
    
    public GTF exclude(GTF input, GTF filter) {   // Sort for speeding up?????
        List t_input = input.getTranscripts();
        List t_filter = filter.getTranscripts();
        List t_output = new ArrayList();
        int offset=0;
        for (int i=0; i<t_input.size(); i++) {
            boolean filter_out=false;
            Transcript t = (Transcript) t_input.get(i);
            while (offset<=t_filter.size()) {
                Transcript f = (Transcript) t_filter.get(offset);
                if (VERBOSE) System.out.print("\n"+t.chr+"("+getChr(t.chr)+"):"+t.start+"-"+t.end+"\t"+f.chr+"("+getChr(f.chr)+"):"+f.start+"-"+f.end+" ("+offset+"/"+t_filter.size()+")");                    
                if ( (getChr(f.chr)>getChr(t.chr)) ||
                     (t.chr.equals(f.chr)&&t.end<f.start)) {
                    break;
                }
                if (t.chr.equals(f.chr) && t.start<=f.end && t.end>=f.start) {
                    filter_out=true;
                    break;
                }
                offset++;
            }
            if (!filter_out) {
                t_output.add(t);
                if (VERBOSE) System.out.print(" *");
            }
        }
        input.setTranscripts(t_output);
        return input;
    }
    
    public void write(GTF output) throws IOException {
        String file=output.getFile().substring(0,output.getFile().lastIndexOf("."))+".output.gtf";
        FileWriter fw = new FileWriter(file);
        BufferedWriter bw = new BufferedWriter(fw);
        List transcripts = output.getTranscripts();
        for (int i=0; i<transcripts.size(); i++) {
            Transcript t = (Transcript) transcripts.get(i);
            bw.write(t.chr+"\t"+t.method+"\t"+t.type+"\t"+t.start+"\t"+t.end+"\t1\t"+t.strand+"\t.\t"+"gene id \""+t.gene+"\"; transcript_id \""+t.id+"\"; exon \""+t.exon+"\"; FPKM \""+t.FPKM+"\"; conf_lo \""+t.conf_lo+"\"; conf_hi \""+t.conf_hi+"\";\n");
        }
        bw.close();
        fw.close();
    }
}
