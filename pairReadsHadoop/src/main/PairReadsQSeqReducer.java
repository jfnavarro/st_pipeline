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
import common.ReadPair;
import common.SequenceId;
import common.WritableMapping;

import org.apache.hadoop.io.Text;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.nio.ByteBuffer;

public class PairReadsQSeqReducer
{
	private static final Log LOG = LogFactory.getLog(PairReadsQSeqReducer.class);

	private Text     outputKey   = new Text();
	private ReadPair outputValue = new ReadPair();

	private ByteBuffer[]      sequence;
	private ByteBuffer[]      quality;
	private WritableMapping[] mapping;

	private int     minBasesThreshold  = 0;
	private boolean dropFailedFilter   = true;
	private boolean warnOnlyIfUnpaired = false;

	private static final byte[] delimByte     = { 9 }; // tab character
	private static final String delim         = "\t";
	private static final char   UnknownBase   = 'N';
	private static final int    INIT_BUF_SIZE = 150;

	public static enum ReadCounters {
		NotEnoughBases,
		FailedFilter,
		Unpaired,
		Dropped
	}

	public void setMinBasesThreshold(int v) { minBasesThreshold = v; }
	public void setDropFailedFilter(boolean v) { dropFailedFilter = v; }
	public void setWarnOnlyIfUnpaired(boolean v) { warnOnlyIfUnpaired = v; }

	public void setup(IMRContext<Text, ReadPair> context)
	{
		sequence = new ByteBuffer[2];
		quality  = new ByteBuffer[2];
		mapping  = new WritableMapping[2];

		for (int i = 0; i < sequence.length; ++i)
		{
			sequence[i] = ByteBuffer.allocate(INIT_BUF_SIZE);
			sequence[i].limit(0).position(0);
			quality[i] = ByteBuffer.allocate(INIT_BUF_SIZE);
			quality[i].limit(0).position(0);
			mapping[i] = new WritableMapping();
		}

		// create counters with a value of 0.
		context.increment(ReadCounters.NotEnoughBases, 0);
		context.increment(ReadCounters.FailedFilter, 0);
		context.increment(ReadCounters.Dropped, 0);
	}

	public void reduce(SequenceId key, Iterable<Text> values, IMRContext<Text,ReadPair> context)
		throws IOException, InterruptedException
	{
		outputKey.set( key.getLocation() );
		outputValue.clear();

		int nReads = 0;
		int nBadReads = 0;
		for (Text read: values)
		{
			++nReads;
			if (nReads > 2)
				throw new RuntimeException("got more than two reads for sequence key " + key + ". Record: " + read);

			int[] fieldsPos = findFields(read);
			// filtered read?
			// If dropFailedFilter is false it shortcuts the test and sets filterPassed directly to true.
			// If it's true then we check whether the field is equal to '1'
			boolean filterPassed = !dropFailedFilter || read.getBytes()[fieldsPos[2]] == (byte)'1';

			if (!filterPassed)
			{
				context.increment(ReadCounters.FailedFilter, 1);
				++nBadReads;
			}
			else if (!checkReadQuality(read, fieldsPos))
			{
				context.increment(ReadCounters.NotEnoughBases, 1);
				++nBadReads;
			}

			// In here we do all the work to prepare the read for output.  It will be written to the
			// appropriate WritableMapping, which will in turn be inserted into the ReadPair outputValue.
			prepMapping(read.getBytes(), fieldsPos, nReads - 1);
		}

		if (nReads == 1)
		{
			context.increment(ReadCounters.Unpaired, nReads);
			if (warnOnlyIfUnpaired)
				LOG.warn("unpaired read!\n" + outputValue.toString());
			else
				throw new RuntimeException("unpaired read for key " + key.toString() + "\nread: " + outputValue.toString());
		}
		else if (nReads != 2)
		{
			throw new RuntimeException("wrong number of reads for key " + key.toString() +
					"(expected 2, got " + nReads + ")\n" + outputValue.toString());
		}

		if (nReads == 2 && nBadReads < nReads) // if they're paired and they're not all bad write. Unpaired are dropped
			context.write(outputKey, outputValue);
		else
			context.increment(ReadCounters.Dropped, nReads);

		context.progress();
	}

	/**
	 * Set a mapping from the mapper-serialized data.
	 *
	 * @param data The byte array from the Text object where the data was serialized by the mapper.
	 * @param fieldPositions Array of positions indexing the byte array, as produced by findFields.
	 * @param index Index of the Mapping we're creating, whether 0 or 1.  We use this to index into
	 * the sequence, quality and mapping arrays and to set the read number in the WritableMapping.
	 *
	 * @post The WritableMapping in mapping[index] will be reset with the data from the input byte
	 * array data.  It will use the ByteBuffers in sequence and quality to store its
	 * sequence and base quality data.  The reset mapping will be set as read 1 or 2 of the
	 * outputValue ReadPair.
	 */
	private void prepMapping(byte[] data, int[] fieldPositions, int index)
	{
		WritableMapping map = mapping[index];
		int length = fieldPositions[1] - 1;
		ensureCapacity(length);

		map.clear();

		// ByteBuffer.clear() resets position to 0 and limit to capacity()
		sequence[index].clear();
		sequence[index].put(data, 0, length).rewind().mark().limit(length);
		map.setSequence(sequence[index]);

		quality[index].clear();
		quality[index].put(data, fieldPositions[1], length).rewind().mark().limit(length);
		map.setBaseQualities(quality[index]);

		if (index == 0)
		{
			map.setIsRead1(true);
			outputValue.setRead1(map);
		}
		else if (index == 1)
		{
			map.setIsRead2(true);
			// set both reads as paired
			mapping[0].setIsPaired(true);
			map.setIsPaired(true);
			outputValue.setRead2(map);
		}
	}

	private void ensureCapacity(int required)
	{
		int current = sequence[0].capacity();
		if (required > current)
		{
			int newSize = required*2;
			// grow
			for (int i = 0; i < sequence.length; ++i)
			{
				ByteBuffer temp;
				// allocate and copy sequence first
				temp = ByteBuffer.allocate(newSize);
				temp.put(sequence[i]).rewind().mark().limit(sequence[i].limit());
				sequence[i] = temp;
				if (mapping[i].getSequence() != null)
					mapping[i].setSequence(sequence[i]);

				// repeat for quality
				temp = ByteBuffer.allocate(newSize);
				temp.put(quality[i]).rewind().mark().limit(quality[i].limit());
				quality[i] = temp;
				if (mapping[i].getBaseQualities() != null)
					mapping[i].setBaseQualities(quality[i]);
			}
		}
	}

	// read format:
	//             read1 <tab> quality1 <tab> filter flag
	// field idx   0           1              2
	private int[] findFields(Text read)
	{
		int[] fieldsPos = new int[3];
		fieldsPos[0] = 0;

		for (int i = 1; i <= 2; ++i)
		{
			fieldsPos[i] = read.find(delim, fieldsPos[i-1]) + 1; // +1 since we get the position of the delimiter
			if (fieldsPos[i] <= 0)
				throw new RuntimeException("invalid read/quality format: " + read.toString());
		}

		int seqLength = fieldsPos[1] - 1;
		int qualLength = fieldsPos[2] - fieldsPos[1] - 1;
		if (seqLength != qualLength)
			throw new RuntimeException("sequence and quality lengths don't match! (got " + seqLength + " and " + qualLength + ")");

		return fieldsPos;
	}

	/**
	 * Verify whether a read satisfies quality standards.
	 * For now this method verifies whether the read has at least
	 * minBasesThreshold known bases (ignoring unknown bases N).
	 */
	protected boolean checkReadQuality(Text read, int[] fieldsPos)
	{
		/* The read's delimiter is at the bytes before the second field starts */
		int readEnd = fieldsPos[1] - 1;

		// The condition is "min number of valid bases".  However, we consider
		// the inverse condition "max number of unknowns".
		// readEnd is also the length of the read fragment
		// readEnd - minBasesThreshold gives us the maximum number of unknowns acceptable.
		int nAcceptableUnknowns = readEnd - minBasesThreshold;

		if (nAcceptableUnknowns < 0) // the fragment is shorter than minBasesThreshold
			return false;

		int nUnknownBases = 0;
		byte[] data = read.getBytes(); // we can work directly in bytes as long as we only have ASCII characters
		for (int pos = 0; pos < readEnd; ++pos)
		{
			if (data[pos] == UnknownBase)
			{
				++nUnknownBases;
				if (nUnknownBases > nAcceptableUnknowns)
					return false;
			}
		}
		return true;
	}
}

