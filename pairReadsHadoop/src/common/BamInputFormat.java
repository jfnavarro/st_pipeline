// Copyright (C) 2011-2012 CRS4.
//
// This file is part of Seal.
//
// Seal is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Seal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with Seal.  If not, see <http://www.gnu.org/licenses/>.

package common;

import common.AbstractTaggedMapping.TagDataType;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.List;

import fi.tkk.ics.hadoop.bam.FileVirtualSplit;
import fi.tkk.ics.hadoop.bam.SAMRecordWritable;
import net.sf.samtools.SAMRecord;

public class BamInputFormat extends FileInputFormat<LongWritable, ReadPair>
{
	private static final Log LOG = LogFactory.getLog(BamInputFormat.class);

	private fi.tkk.ics.hadoop.bam.BAMInputFormat bamImpl;

	public static class BamRecordReader extends RecordReader<LongWritable, ReadPair>
	{
		private ReadPair value;
		private WritableMapping mapping;
		private RecordReader<LongWritable, SAMRecordWritable> rrImpl;

		public BamRecordReader(fi.tkk.ics.hadoop.bam.BAMRecordReader finBamRR)
		{
			rrImpl = finBamRR;
		}

		@Override
		public void initialize(InputSplit genericSplit, TaskAttemptContext context) throws IOException, InterruptedException
		{
			rrImpl.initialize(genericSplit, context);

			value = new ReadPair();
			mapping = new WritableMapping();
		}

		@Override
		public void close() throws IOException
		{
			rrImpl.close();
		}

		@Override
		public float getProgress() throws IOException, InterruptedException { return rrImpl.getProgress(); }

		@Override
		public LongWritable getCurrentKey() throws IOException, InterruptedException { return rrImpl.getCurrentKey(); }

		@Override
		public ReadPair getCurrentValue() { return value; }

		@Override
		public boolean nextKeyValue() throws IOException, InterruptedException
		{
			if (rrImpl.nextKeyValue())
			{
				value.clear();
				SAMRecord sam = rrImpl.getCurrentValue().get();
				// copy data from SAMRecord to our mapping
				readSamRecord(sam, mapping);
				if (mapping.isRead2())
					value.setRead2(mapping);
				else // anything that's not explicitly labelled as "read 2" goes in as read 1.
					value.setRead1(mapping);
				return true;
			}
			else
				return false;
		}

		protected void readSamRecord(SAMRecord sam, WritableMapping mapping)
		{
			mapping.clear();

			mapping.setName(sam.getReadName());
			mapping.setSequence(ByteBuffer.wrap(sam.getReadBases()));
			mapping.setFlag(sam.getFlags());

			// we have to map base qualities from raw phred-scaled scores to
			// Phred+33, as per our convention
			byte[] originalQ = sam.getBaseQualities();
			ByteBuffer q = ByteBuffer.allocate(originalQ.length);
			byte[] newQ = q.array();
			int newQStart = q.position();
			for (int i = 0; i < originalQ.length; ++i)
			{
				if (originalQ[i] < 0 || originalQ[i] > Utils.SANGER_MAX)
				{
					throw new FormatException(
					  "base quality score out of range for BAM format (found " + originalQ[i] +
					  " but acceptable range is [0," + Utils.SANGER_MAX + "]).");
				}
				newQ[ newQStart + i ] = (byte)(originalQ[i] + Utils.SANGER_OFFSET);
			}
			q.rewind().mark();
			mapping.setBaseQualities(q);

			// now get the rest of the fields
			if (mapping.isMapped())
			{
				mapping.setContig(sam.getReferenceName());
				mapping.set5Position(sam.getAlignmentStart());
				mapping.setMapQ(sam.getMappingQuality());
				mapping.setAlignment(AlignOp.scanCigar(sam.getCigarString()));
			}

			if (mapping.isPaired() && mapping.isMateMapped())
				mapping.setTemplateLength(Math.abs(sam.getInferredInsertSize()));

			for (SAMRecord.SAMTagAndValue tag: sam.getAttributes())
			{
				Object value = tag.value;
				TagDataType type = null;

				if (value instanceof String)
					type = TagDataType.String;
				else if (value instanceof Character)
					type = TagDataType.Char;
				else if (value instanceof Float)
					type = TagDataType.Float;
				else if (value instanceof Integer || value instanceof Byte || value instanceof Short)
					type = TagDataType.Int;
				else
				{
					LOG.warn("dropping tag of unsupported type " + value.getClass().toString());
					continue;
				}

				mapping.setTag(tag.tag, type, value.toString());
			}
		}
	}

	public BamInputFormat()
	{
		bamImpl = new fi.tkk.ics.hadoop.bam.BAMInputFormat();
	}

	@Override
	public List<InputSplit> getSplits(JobContext job) throws IOException
	{
		return bamImpl.getSplits(job);
	}

	/**
	 * Public method to create FileVirtualSplits.
	 *
	 * Useful for unit testing.
	 */
	public List<InputSplit> getVirtualSplits(List<InputSplit> fileSplits, Configuration conf) throws IOException
	{
		return bamImpl.getSplits(fileSplits, conf);
	}

	@Override
	public RecordReader<LongWritable,ReadPair> createRecordReader(InputSplit split, TaskAttemptContext context)
		throws IOException, InterruptedException
	{
		return new BamRecordReader((fi.tkk.ics.hadoop.bam.BAMRecordReader)bamImpl.createRecordReader(split, context));
	}

	@Override
	protected boolean isSplitable(JobContext context, Path file)
	{
		return bamImpl.isSplitable(context, file);
	}
}
