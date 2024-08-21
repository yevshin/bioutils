package yevshin;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BWAAlignStats {
	
	static int mapq = 10;
	
	static class PairAlign
	{
		SAMRecord r1, r2;
		List<SAMRecord> suppR1, suppR2;
		
		void add(SAMRecord r)
		{
			checkExpectations(r);
			
			//BWA outputs chimeric alignments as supplementary (when read is split into parts)
			if(r.getSupplementaryAlignmentFlag())
				addSupplementary(r);
			else
			{
				if(r.getFirstOfPairFlag())
				{
					if(r1 != null)
						throw new RuntimeException("Duplicated: /// " + r + " /// " + r1);
					r1 = r;
				}
				else
				{
					if(r2 != null)
						throw new RuntimeException("Duplicated: /// " + r + " /// " + r2);
					r2 = r;
				}
			}
			
		}

		private void addSupplementary(SAMRecord r) {
			if(r.getFirstOfPairFlag())
			{
				if( suppR1 == null )
					suppR1 = new ArrayList<>();
				suppR1.add(r);
			}
			else
			{
				if( suppR2 == null )
					suppR2 = new ArrayList<>();
				suppR2.add(r2);
			}
		}

		private void checkExpectations(SAMRecord r) {
			if(!r.getReadPairedFlag())
				throw new RuntimeException("Not paired: " + r.getSAMString());
			if(r.getNotPrimaryAlignmentFlag())
				throw new RuntimeException("Not primary: " + r.getSAMString());
			if(r.getReadFailsVendorQualityCheckFlag())
				throw new RuntimeException("Vendor quality: " + r.getSAMString());
			if(r.getDuplicateReadFlag())
				throw new RuntimeException("PCR duplicate: " + r.getSAMString());
		}
	}
	
	public static void main(String[] args) throws IOException {
		String inbam = args[0];//should be direct output from BWA-MEM paired-end, default params
		if(args.length > 1)
			mapq = Integer.parseInt(args[1]);
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamInputResource in = SamInputResource.of(new File(inbam));
		SamReader reader = factory.open(in);
		
		SAMRecordIterator it = reader.iterator();
		if(!it.hasNext())
			return;
		SAMRecord r = it.next();
		
		OutputStream novelOutR1 = new GZIPOutputStream(new FileOutputStream("novel_sequences_R1.fq.gz"));
		OutputStream novelOutR2 = new GZIPOutputStream(new FileOutputStream("novel_sequences_R2.fq.gz"));
		novelSequencesR1 = new BufferedWriter(new OutputStreamWriter(novelOutR1));
		novelSequencesR2 = new BufferedWriter(new OutputStreamWriter(novelOutR2));
	
		boolean nextPair = true;
		while(nextPair)
		{
			PairAlign p = new PairAlign();
			p.add(r);
			String name = r.getReadName();
			
			nextPair = false;
			while(it.hasNext())
			{
				r = it.next();
				if(r.getReadName().equals(name))
					p.add(r);
				else
				{
					nextPair = true;
					break;
				}
			}
			
			processPair(p);
		}
		
		reader.close();
		
		novelSequencesR1.close();
		novelSequencesR2.close();
		System.out.println("Total pairs: " + S.total);
		System.out.println("Mapped uniquely: " + S.mapped + " (" + String.format("%.2f", S.mapped*100.0/S.total) + "%)" );
		System.out.println("Chimeric (structural variants or PCR artifacts): " + S.chimeric + " (" + String.format("%.2f", S.chimeric*100.0/S.total) + "%)" );
		System.out.println("Ambigous mapping due to repeats: " + S.ambiguosMapping + " (" + String.format("%.2f", S.ambiguosMapping*100.0/S.total) + "%)" );
		System.out.println("Unmapped due to low base call quality: " + S.lowBaseCallQuality + " (" + String.format("%.2f", S.lowBaseCallQuality*100.0/S.total) + "%)" );
		System.out.println("Novel sequences, partially mapped: " + S.novelSequencePartiallyMapped + " (" + String.format("%.2f", S.novelSequencePartiallyMapped*100.0/S.total) + "%)" );
		System.out.println("Novel sequences, fully unmapped: " + S.novelSequenceFullyUnmapped + " (" + String.format("%.2f", S.novelSequenceFullyUnmapped*100.0/S.total) + "%)" );
	}
	
	static class Stat
	{
		int total;
		int mapped;
		
		//unmapped
		int lowBaseCallQuality;//due to sequencer call quality
		int ambiguosMapping;//due to repeats
		
		//Sequence absent in reference
		int novelSequencePartiallyMapped;
		int novelSequenceFullyUnmapped;
		
		int chimeric;//due to structural genomic rearrangments or PCR artifacts
	}
	
	static Stat S = new Stat();

	private static void processPair(PairAlign p) throws IOException {
		S.total++;
		
		if(p.r1 == null || p.r2 == null)
			throw new RuntimeException("Missing r1, r2");
		
		if(p.r1.getReadUnmappedFlag() || p.r2.getReadUnmappedFlag())
		{
			if(lowBaseCallQual(p.r1) || lowBaseCallQual(p.r2))
				S.lowBaseCallQuality++;
			else
			{
				if(p.r1.getReadUnmappedFlag() && p.r2.getReadUnmappedFlag())
				{
					S.novelSequenceFullyUnmapped++;
					reportNovelSequences(p.r1, p.r2);
				}
				else
					S.novelSequencePartiallyMapped++;
			}
		}else
		{
			//both mapped 
			if(p.r1.getMappingQuality() < mapq || p.r2.getMappingQuality() < mapq)//check tags XA, SA
				S.ambiguosMapping++;
			else
			{
				if(p.suppR1 != null || p.suppR2 != null)
					S.chimeric++;
				else
				{
					if(fractionUnmapped(p.r1) >= 0.2 || fractionUnmapped(p.r2) >= 0.2)
						S.novelSequencePartiallyMapped++;
					else
					{
						if(!p.r1.getProperPairFlag() || !p.r2.getProperPairFlag() || p.r1.getReferenceIndex() != p.r2.getReferenceIndex() )
							S.chimeric++;
						else
							S.mapped++;
					}
				}
			}
		}
		
	}

	static BufferedWriter novelSequencesR1;
	static BufferedWriter novelSequencesR2;
	private static void reportNovelSequences(SAMRecord r1, SAMRecord r2) throws IOException {
		printFastqRecord(novelSequencesR1, r1);
		printFastqRecord(novelSequencesR2, r2);
	}

	public static void printFastqRecord(BufferedWriter writer, SAMRecord r1) throws IOException {
		String seq = r1.getReadString();
		String qual = r1.getBaseQualityString();
		if(r1.getReadNegativeStrandFlag())
		{
			seq = rc(seq);
			qual = new StringBuilder(qual).reverse().toString();
		}
		writer
		.append('@').append( r1.getReadName() ).append('\n')
		.append(seq).append('\n');
		
		writer.append("+\n").append(qual).append('\n');
	}
	
	
	public static String rc(String seq)
    {
        char[] res = new char[seq.length()];
        for(int i = 0; i < seq.length(); i++)
        {
            char c = seq.charAt( i );
            switch(c)
            {
                case 'A': res[res.length - i - 1] = 'T'; break;
                case 'C': res[res.length - i - 1] = 'G'; break;
                case 'G': res[res.length - i - 1] = 'C'; break;
                case 'T': res[res.length - i - 1] = 'A'; break;
                case 'N': res[res.length - i - 1] = 'N'; break;
                default:
                    throw new AssertionError();
            }
        }
        return new String(res);
    }

	private static double fractionUnmapped(SAMRecord r) {
		int n = 0;
		for(CigarElement e : r.getCigar())
		{
			CigarOperator op = e.getOperator();
			if(op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.I)
				n += e.getLength();
		}
		return ((double)n) / r.getReadLength();
	}

	private static boolean lowBaseCallQual(SAMRecord r) {
		byte[] qs =  r.getBaseQualities();
		int low = 0;
		for(int i = 0; i < qs.length; i++)
		{
			if(qs[i] < 20)
				low++;
		}
		
		return low > qs.length/2;
	}
}
