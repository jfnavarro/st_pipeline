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

package main;

import common.IMRContext;
import common.SequenceId;

import fi.tkk.ics.hadoop.bam.SequencedFragment;

import org.apache.hadoop.io.Text;

import java.io.IOException;

public class PairReadsQSeqMapper
{
	private StringBuilder builder;
	private SequenceId sequenceKey = new SequenceId();
	private Text sequenceValue = new Text();

	private static final int LINE_SIZE = 500;
	private static final byte[] Delim = { 9 }; // tab
	private static final byte[] ZeroOne = { '0', '1' };

	private boolean makeTraditionalIds = false;

	public void setMakeTraditionalIds(boolean v) {
		makeTraditionalIds = v;
	}

	public void setup()
	{
		builder = new StringBuilder(LINE_SIZE);
	}

	public void map(Text readId, SequencedFragment read, IMRContext<SequenceId, Text> context) throws IOException, InterruptedException
	{
		// build the key
		builder.delete(0, builder.length());

		// field up and including the index number goes in the location.  The read is on its own.
		if (read.getRead() == null)
			throw new RuntimeException("Cannot get read number from read: " + readId);

		if (read.getLane() != null && read.getTile() != null && read.getXpos() != null && read.getYpos() != null)
		{
			appendIdToBuilder(builder, read); // appends the read id to the builder provided
			// finally the index field
			builder.append("#").append(read.getIndexSequence() == null ? '0' : read.getIndexSequence());
			sequenceKey.set(builder.toString(), read.getRead());
		}
		else
		{
			// maybe it's a fastq id with a trailing read number (/1 or /2)
			if (readId.getLength() > 2)
			{
				int last = readId.getLength() - 1;
				if (readId.charAt(last - 1) == '/')
				{
					// truncate the /[12] from the read id
					// last == length - 1.  We want length - 2 bytes, which is equal to last - 1
					sequenceKey.set(Text.decode(readId.getBytes(), 0, last - 1), read.getRead());
				}
				else
					throw new RuntimeException("Didn't find /read_number at end of the read id.  Please use qseq files or fastq with illumina-formatted name tags.");
			}
			else
				throw new RuntimeException("Read id " + readId + " is too short.   Please use qseq files or fastq with illumina-formatted name tags.");
		}

		// then the tab-delimited value
		sequenceValue.clear();
		sequenceValue.append(read.getSequence().getBytes(), 0, read.getSequence().getLength());
		sequenceValue.append(Delim, 0, Delim.length);
		sequenceValue.append(read.getQuality().getBytes(), 0, read.getQuality().getLength());
		sequenceValue.append(Delim, 0, Delim.length);
		// the filter flag is optional.  If it's absent we assume the read passes filtering.
		sequenceValue.append(ZeroOne, (read.getFilterPassed() == null || read.getFilterPassed() ? 1 : 0), 1);

		context.write(sequenceKey, sequenceValue);
		context.progress();
	}

	protected StringBuilder appendIdToBuilder(StringBuilder builder, SequencedFragment read)
	{
		builder.append(read.getInstrument() == null ? "" : read.getInstrument());
		if (makeTraditionalIds)
		{
			builder.append("_").append(read.getRunNumber() == null ? "" : read.getRunNumber());
		}
		else
		{
			builder.append(":").append(read.getRunNumber() == null ? "" : read.getRunNumber());
			builder.append(":").append(read.getFlowcellId() == null ? "" : read.getFlowcellId());
		}

		builder.append(":").append(read.getLane());
		builder.append(":").append(read.getTile());
		builder.append(":").append(read.getXpos());
		builder.append(":").append(read.getYpos());

		return builder;
	}
}
