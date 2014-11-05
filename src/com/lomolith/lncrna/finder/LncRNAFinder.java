/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.lncrna.finder;

import com.lomolith.annotation.GTFManager;
import com.lomolith.common.model.GTF;
import com.lomolith.common.model.Transcript;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author mittens
 */
public class LncRNAFinder {

    private String file="transcripts.gtf";
    private String remove="";
    private String annot="genes.gtf";
    private String reference="hg19.fa";
    private String dir="";
    
    private boolean paired=true;
    private int length=300;
    
    public boolean VERBOSE=false;   // Not implemented yet
    public boolean WRITE=false;      // Not implemented yet
    
    public LncRNAFinder() {
    }

    public LncRNAFinder(String workDir) {
        dir=workDir;
    }
    
    public List getTrainSequences(String refGTF, String refFA, GTF annotation) throws IOException {
        reference=refFA;
        annot=refGTF;
        return getTrainSequences(false, annotation);
    }
    
    public List getTrainSequences(String refGTF, String refFA, boolean filter, GTF annotation) throws IOException {
        reference=refFA;
        annot=refGTF;
        return getTrainSequences(filter, annotation);
    }
    
    public List getTrainSequences(String refGTF, String refFA) throws IOException {
        reference=refFA;
        annot=refGTF;
        return getTrainSequences();
    }
    
    public List getTrainSequences(String refGTF) throws IOException {
        annot=refGTF;
        return getTrainSequences(true);
    }
    
    public List getTrainSequences(String refGTF, boolean filter) throws IOException {
        annot=refGTF;
        return getTrainSequences(filter);
    }
    
    public List getTrainSequences() throws IOException {
        return getTrainSequences(false);
    }
    
    public List getTrainSequences(boolean sizeFilter) throws IOException {
        return getTrainSequences(sizeFilter, null);
    }
    
    public List getTrainSequences(boolean sizeFilter, GTF exist_annot) throws IOException {
        GTFManager gm = new GTFManager();
        gm.setFile(dir+"/"+annot);
        gm.FILTER_SIZE=sizeFilter;
        System.out.print("Reading annotation: "+annot+"...");
        GTF annotation = gm.read(exist_annot);
        System.out.println("done. Total "+annotation.getSize()+" transcripts.");
        System.out.print("Retrieving sequences for "+annotation.getSize()+" transcripts...");
        annotation.reference=dir+"/"+reference;
        annotation.getSeq();
        System.out.println("done.");
        if (WRITE) gm.write(annotation);
        return annotation.getTranscripts();
    }
    
    public List getTriFrequency() throws IOException {
        GTFManager gm = new GTFManager();
        gm.setFile(dir+"/"+annot);
        System.out.print("Reading annotation: "+annot+"...");
        GTF annotation = gm.read();
        if (WRITE) gm.write(annotation);
        System.out.println("done. Total "+annotation.getSize()+" transcripts.");

        gm.setFile(dir+"/"+file);
        System.out.print("Reading target : "+file+"...");
        GTF gtf = gm.read();
        if (WRITE) gm.write(gtf);
        System.out.println("done. Total "+gtf.getSize()+" transcripts.");

        System.out.print("Removing annotated regions...");
        gtf = gm.exclude(gtf, annotation);
        System.out.println("done. Total "+gtf.getSize()+" transcripts remain.");

        System.out.print("Retrieving sequences for "+gtf.getSize()+" transcripts...");
        gtf.reference=dir+"/"+reference;
        gtf.getSeq();
        System.out.println("done.");

        int start_codon=0;
        int orf=0;
        List detected = gtf.getTranscripts();
        List candidates = new ArrayList();

        for (int i=0; i<detected.size(); i++) {
            Transcript t = (Transcript) detected.get(i);
            if ( t.sequence.toUpperCase().indexOf("ATG")!=-1) {
//                    System.out.printf("Found start codon: %3d %2s:%10d-%10d %6d %6d %10s\n",i,t.chr,t.start,t.end,t.end-t.start+1,t.sequence.lengthORF(), t.sequence.substring(0,10));
                start_codon++;
                String sub = t.sequence.substring(t.sequence.toUpperCase().indexOf("ATG"));

                int pos=findStopCodon(sub);
                if (pos!=-1 && pos+3>length) {
                    orf++;
                    t.orf_seq=sub.substring(0,pos+3);
                    t.orf_open="P";
                    detected.set(i, t);
                }
                else {
                    if (!paired) {
                        orf++;
                        t.orf_seq=sub;
                        t.orf_open="S";
                        detected.set(i, t);
                    }
                }
            }
            else {
                int ochre=t.sequence.toUpperCase().indexOf("TAA");
                int amber=t.sequence.toUpperCase().indexOf("TAG");
                int opal_umber=t.sequence.toUpperCase().indexOf("TGA");
                int pos=Math.min(Math.min(ochre, amber), opal_umber);
                if ( pos!=-1 && pos+3>length && !paired) {
                    orf++;
                    t.orf_seq=t.sequence.substring(0,pos+3);
                    t.orf_open="E";
                    detected.set(i, t);
                };
            }
        }

        for (Object o: detected) {
            Transcript t = (Transcript) o;
            if (t.orf_open==null) {
                System.out.println(t.sequence.length());
                System.out.println(t.sequence);
                candidates.add(t);
            }
        }
        System.out.println("Found orf in "+orf+" transcripts from total "+detected.size());
        System.out.println("Found start codon only in "+(start_codon-orf)+" transcripts from total "+detected.size());
        System.out.println("Found "+candidates.size()+" candidates from total "+detected.size());
        
        return candidates;
    }
    
    private int findStopCodon(String seq) {
        for (int i=0; i<seq.length()-3; i=i+3) {
            String sub=seq.substring(i,i+3).toUpperCase();
            if (sub.equals("TAA")||sub.equals("TAG")||sub.equals("TGA")) return i;
        }
        return -1;
    }
}
