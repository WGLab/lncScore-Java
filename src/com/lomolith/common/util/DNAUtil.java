/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.common.util;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;

public class DNAUtil {

    public boolean CONNECTION=true;
    private String strUrl = "http://www.ensembl.org/das/Homo_sapiens.GRCh37.reference/sequence?segment=";

    public DNAUtil() {
        try {
            URL url=new URL("http://www.google.com");
            URLConnection con=url.openConnection();
            con.getInputStream();    
        }
        catch (Exception e) {
            CONNECTION=false;
            System.err.println("\nWarning: Internet disconnected.");
        }
    }
    
    public String getSeq(String coord) throws MalformedURLException, IOException {
        if (CONNECTION && coord.indexOf(":")!=-1 && coord.indexOf("-")!=-1) {
            String chr=coord.substring(0,coord.indexOf(":"));
            long start=Long.parseLong(coord.substring(coord.indexOf(":")+1,coord.indexOf("-")));
            long end=Long.parseLong(coord.substring(coord.indexOf("-")+1));
            return getSeq(chr, start, end);
        }
        else return null;
    }
    
    public String getSeq(String chr, long start, long end) throws MalformedURLException, IOException {
        String seq=null;
        if (CONNECTION) {
            strUrl=strUrl+chr.toUpperCase().replaceAll("CHR","")+":"+start+","+end;
            BufferedInputStream in;
            StringBuffer sb=new StringBuffer();
            URL url = new URL(strUrl);
            URLConnection urlConnection = url.openConnection();
            in = new BufferedInputStream(urlConnection.getInputStream());

            byte[] bufRead = new byte[4096];
            int lenRead = 0;
            while ((lenRead = in.read(bufRead)) > 0)
                sb.append(new String(bufRead, 0, lenRead));
            seq=sb.toString().replaceAll("<[^>]*>", "").trim();
        }
        return seq;
    }
}
