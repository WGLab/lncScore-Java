package com.lomolith.common.model;

import com.lomolith.sequence.FastaReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author mittens
 */
public class GTF {
    private String filename="";
    public List transcripts=new ArrayList();

    public boolean FILTER_SIZE=true;
    public String reference="/home/jihongkim/workspace/tade/hg19.fa";
    
    public GTF(String name) {
        filename=name;
    }
    
    public void setFile(String name) {
        filename=name;
    }
    public String getFile() {
        return filename;
    }
    public int getSize() {
        return transcripts.size();
    }
    
    public List getTranscript(String name) { 
        List returnValue = new ArrayList();
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (t.gene.equals(name) || t.id.equals(name)) returnValue.add(t);
        }
        return returnValue;
    }
    public List getTranscripts(String chr, double start, double end) {
        List returnValue = new ArrayList();
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (t.chr.equals(chr) && t.start<=end && t.end>=start) returnValue.add(t);
        }
        return returnValue;
    }
    public List getTranscripts(float fpkm) {
        List returnValue = new ArrayList();
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (t.FPKM>=fpkm) returnValue.add(t);
        }
        return returnValue;
    }
    public List getTranscripts(String type) {
        List returnValue = new ArrayList();
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (t.type.equals(type)) returnValue.add(t);
        }
        return returnValue;
    }
    
    public void setTranscript(String line) {                                    // Currently, it ignores duplicates of transcript ID.
        String[] cols=line.split("\t");
        if (cols[2].equals("gene")||cols[2].equals("transcript")) {
            Transcript t = new Transcript();
            t.chr=cols[0].toUpperCase().replaceAll("CHR","");
            t.start=Long.parseLong(cols[3]);
            t.end=Long.parseLong(cols[4]);
            if (t.chr.indexOf("_")==-1) {
                if ((t.end-t.start+1>200 && FILTER_SIZE) || !FILTER_SIZE) {
                    t.method=cols[1];
                    t.type=cols[2];
                    t.strand=cols[6].trim();
                    String[] desc=cols[8].split(";");
                    for (int i=0; i<desc.length; i++) {
                        String v = desc[i].trim().substring(desc[i].trim().indexOf(" ")+1).replaceAll("\"","");
                        if (desc[i].contains("gene_id")) t.gene=v;
                        else if (desc[i].contains("transcript_id")) t.id=v;
                        else if (desc[i].contains("FPKM")) t.FPKM=Float.parseFloat(v);
                        else if (desc[i].contains("conf_hi")) t.conf_hi=Float.parseFloat(v);
                        else if (desc[i].contains("conf_lo")) t.conf_lo=Float.parseFloat(v);
                        else if (desc[i].contains("exon_number")) t.exon=Integer.parseInt(v);
                    }
                    transcripts.add(t);
                }
            }
        }
        else if (cols[2].equals("exon")) {                                      // It assume every exons should be appeared after apearance of transcripts.
            long start=Long.parseLong(cols[3]);
            long end=Long.parseLong(cols[4]);
            String[] desc=cols[8].split(";");
            for (int i=0; i<desc.length; i++) {
                String v = desc[i].trim().substring(desc[i].trim().indexOf(" ")+1).replaceAll("\"","");
                if (desc[i].contains("transcript_id")) {
                    for (int c=0; c<transcripts.size(); c++) {
                        Transcript t = (Transcript) transcripts.get(c);
                        if (t.id.equals(v)) {
                            t.size+=end-start+1;
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    
    public List getTranscripts() {
        return transcripts;
    }
    
    public void sortTranscripts() {
        sortTranscripts(true);
    }
    
    public void sortTranscripts(final boolean asc) {
        Collections.sort(transcripts, new Comparator(){
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
            public int compare(Object left, Object right){
                Transcript leftProbe = (Transcript)left;
                Transcript rightProbe = (Transcript)right;

                double leftValue = leftProbe.start;
                double rightValue = rightProbe.start;
                double leftValue2 = leftProbe.end;
                double rightValue2 = rightProbe.end;
                int leftChr = getChr(leftProbe.chr);
                int rightChr = getChr(rightProbe.chr);
                if (asc) {
                    if (!leftProbe.chr.equals(rightProbe.chr)) {
                        if (leftChr<rightChr) return -1;
                        else return 1;
                    }
                    else if (leftValue<rightValue) return -1;
                    else if (leftValue>rightValue) return 1;
                    else if (leftValue2<rightValue2) return -1;
                    else if (leftValue2>rightValue2) return 1;
                    else return 0;
                }
                else {
                    if (leftChr!=rightChr) {
                        if (leftChr<rightChr) return 1;
                        else return -1;
                    }
                    else if (leftValue>rightValue) return -1;
                    else if (leftValue<rightValue) return 1;
                    else if (leftValue2>rightValue2) return -1;
                    else if (leftValue2<rightValue2) return 1;
                    else return 0;
                }
            }
        });
    }
        
    public void setTranscripts(List t) {
        transcripts = new ArrayList(t);
        sortTranscripts();
    }
    
    public void removeTranscript(Transcript target) {
        List remains = new ArrayList();
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (!t.chr.equals(target.chr) || t.start>target.end || t.end<target.start) remains.add(t);
        }
        transcripts = new ArrayList(remains);
    }
    
    public boolean findTranscript(Transcript target) {
        //Object z[] = transcripts.keySet().toArray();
        //for (int i=0; i<z.length; i++) System.out.println(z[i]+" / ");
        for (Object k:transcripts) {
            Transcript t = (Transcript) k;
            if (t.chr.equals(target.chr) && t.start<=target.end && t.end>=target.start) return true;
        }
        return false;
    }
    
    public void getSeq() throws IOException {
        FastaReader f=new FastaReader(new File(reference));
//        DNAUtil du = new DNAUtil();
        for (int i=0; i<transcripts.size(); i++) {
            Transcript t = (Transcript) transcripts.get(i);
//            if (du.CONNECTION) t.sequence = du.getSeq(t.chr, t.start, t.end);
            byte[] temp_seq = f.getSequence("chr"+t.chr, (int)t.start, (int)(t.end-t.start+1));
            if (temp_seq!=null) {
                t.sequence = new String(temp_seq);
                if (t.strand.equals("-")) t.sequence = new StringBuffer(t.sequence).reverse().toString();
            }
        }
        List filtered = new ArrayList();
        for (int i=0; i<transcripts.size();i++) {
            Transcript t = (Transcript) transcripts.get(i);
            if (t.sequence!=null && !t.sequence.equals("")) filtered.add(t);
        }
        transcripts = filtered;
    }
}
