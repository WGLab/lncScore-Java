package com.lomolith.sequence;

import com.lomolith.common.model.FastaIndex;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.Map;

public class FastaReader {
    private Map<String, FastaIndex> name2index=new HashMap<String, FastaIndex>();
    private RandomAccessFile io=null;

    public FastaReader(File fastaFile) throws IOException {
        File fastaIndex=new File(fastaFile.getPath()+".fai");
        BufferedReader in=new BufferedReader(new FileReader(fastaIndex));
        String line;
        while((line=in.readLine())!=null){
            String tokens[]=line.split("\t");
            FastaIndex index=new FastaIndex();
            index.names.add(tokens[0]);
            index.source=fastaFile;
            index.sequence_length=Integer.parseInt(tokens[1]);
            index.seqStart=Long.parseLong(tokens[2]);
            index.lineSize=Integer.parseInt(tokens[3]);
            this.name2index.put(tokens[0],index);
        }
        in.close();
    }

    public byte[] getSequence(String chromosome, int start, int length) throws IOException {
        try {
            FastaIndex index=name2index.get(chromosome);
            if(index==null) throw new IOException("Sequence \""+chromosome+"\" was not indexed");
            if(io==null) io=new RandomAccessFile(index.getFile(), "r");
            return index.getSequence(io, start, length);
        }
        catch (IllegalArgumentException e) {
            return null;
        }
    }

    public void close() throws IOException {
        try { 
            if(io!=null) io.close(); 
        } 
        catch(IOException err) {
            err.printStackTrace();
        }
        io=null;
    }

    public static void main(String[] args) {
        try {
            FastaReader f=new FastaReader(new File(System.getenv("HG18")));
            System.err.println(new String(f.getSequence("chrM", 0, 100)));
            System.err.println(new String(f.getSequence("chrM", 2, 100)));
        } 
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}