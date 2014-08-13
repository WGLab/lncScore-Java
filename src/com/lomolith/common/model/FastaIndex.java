package com.lomolith.common.model;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashSet;
import java.util.Set;

/**
 * Source by   https://code.google.com/u/105051423568778287795/
 *        from https://code.google.com/p/code915/
 * 
 * FastaIndex describes a sequence in a FASTA file. what is its name, where does the sequence starts, what
 * is the length of the lines. It allows to quickly read a sub-sequence in a fasta file.
 * @see IndexedGenome
 * @see FastaIndexer
 */

public class FastaIndex
        {
        private static int BUFFER_SIZE=5096;
        /** FASTA file as source */
        public File source;
        
        protected long nameStart=-1L;
        protected long nameEnd=-1L;
        public long seqStart=-1L;
        protected long seqEnd=-1L;
        /** size of the lines in the FASTA file*/
        public int lineSize=0;
        /** name(s) for this sequence */
        public Set<String> names=new HashSet<String>();
        /** size of this sequence */
        public int sequence_length=0;
        
        
        /** fasta file for this sequence */
        public File getFile()
                {
                return source;
                }
        /** name(s) for this sequence e.g. 'chr1, chr01, 01, 1' .... */
        public Set<String> getNames()
                {
                return names;
                }
        /** size of this sequence */
        public int getLength()
                {
                return sequence_length;
                }
        
        /** get file offset for the first base of this sequence */
        public long getSeqStart() {
                return seqStart;
        }
        
        /** get file offset after the last base of this sequence */
        public long getSeqEnd()
                {
                return seqEnd;
                }
        /** set the size of the BUFFER when involing  'getSequence' */
        public static void setBufferSize(int buffSize)
                {
                if(buffSize<=0) throw new IllegalArgumentException("buffSize<=0 : "+buffSize);
                BUFFER_SIZE = buffSize;
                }
        
        /** get the size of the BUFFER when involing  'getSequence' */
        public static int getBufferSize()
                {
                return BUFFER_SIZE;
                }
        
        /** 
         * returns a sequence from a RandomAccessFile for this index
         * @param in a fasta file opened with a RandomAccessFile
         * @param start the start offset of the sequence (first base=0)
         * @param length number of symbols to read
         * @return an array of bytes containing the sequence
         * @throws IOException
         */
        public byte[] getSequence(
                                RandomAccessFile in,
                                int start,
                                int length
                                ) throws IOException
                {
                if(start<0) throw new IllegalArgumentException("start<0:"+start);
                if(start+length> getLength()) throw new IllegalArgumentException("start="+start+" length="+length+"> size="+getLength());
                /* the array of byte where we store the final sequence */
                byte seq[]=new byte[length];
                
                /** the fasta row index */
                int row_index=start/this.lineSize;
                /** index of 'start' in this current row */
                int index_in_row=start-row_index*this.lineSize;
                /** prepare a buffer to read some bytes from the fasta file */
                byte  buffer[]=new byte[FastaIndex.getBufferSize()];
                /**  number of bytes in the buffer */
                int buffer_length=0;
                /** current position in the buffer */
                int index_in_buffer=0;
                /** current position in the final array of bytes (sequence) */
                int index_in_seq=0;
                /** move the IO cursor to the beginning of the sequence */
                in.seek(
                                this.getSeqStart()+
                                row_index*(this.lineSize+1)+
                                index_in_row
                                );
                while(length>0)
                        {
                        /* buffer empty ? fill it */
                        if(buffer_length==0)
                                {
                                buffer_length=in.read(buffer);
                                index_in_buffer=0;
                                if(buffer_length<=0) throw new IOException("cannot fill buffer");
                                }
                        /* number of byte to copy into the sequence */
                        int n_to_copy= Math.min(buffer_length-index_in_buffer,Math.min(length,this.lineSize-index_in_row));
                        //System.err.println("n="+buffer.length+" L="+buffer_length+"\t"+index_in_buffer+" \t"+n_to_copy);
                        //System.err.println("n_to_copy="+n_to_copy+" index_in_row="+index_in_row+" \""+new String(buffer,index_in_buffer,n_to_copy)+"\"");
                        
                        /* copy the bytes from the buffer to the sequence */
                        System.arraycopy(buffer, index_in_buffer, seq, index_in_seq, n_to_copy);
                        
                        /* move the offsets */
                        index_in_seq+=n_to_copy;
                        length-=n_to_copy;
                        index_in_buffer+=n_to_copy;
                        index_in_row=(index_in_row+n_to_copy)%this.lineSize;
                        
                        /* check the next  input is a CR/LF */  
                        if(length>0)
                                {
                                if(index_in_row==0)
                                        {
                                        /* buffer is filled, read from input stream */
                                        if(index_in_buffer==buffer_length)
                                                {
                                                if(in.read()!='\n')
                                                        {
                                                        throw new IOException("I/O error expected a carriage return at offset "+in.getFilePointer());
                                                        }
                                                }
                                        else /* check the buffer */
                                                {
                                                if(buffer[index_in_buffer]!='\n')
                                                        {
                                                        throw new IOException("I/O error expected a carriage return in buffer at "+index_in_buffer);
                                                        }
                                                index_in_buffer++;
                                                }
                                        }
                                        
                                /* reset the buffer it is filled */
                                if(index_in_buffer==buffer_length)
                                        {
                                        buffer_length=0;
                                        index_in_buffer=0;
                                        }
                                }
                        
                        
                        }
                
                

                return seq;
                }
        
        
        
        }