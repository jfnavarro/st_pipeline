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

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.mapreduce.lib.input.LineRecordReader;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;

import java.io.IOException;

public class SamInputFormat extends FileInputFormat<LongWritable, ReadPair>
{
	public static class SamRecordReader extends RecordReader<LongWritable, ReadPair>
	{
		private LineRecordReader lineReader;
		private ReadPair value;
		private FileSplit split; // memorize it for error messages

		public void initialize(InputSplit genericSplit, TaskAttemptContext context) throws IOException
		{
			lineReader = new LineRecordReader();
			lineReader.initialize(genericSplit, context);

			split = (FileSplit)genericSplit;

			value = new ReadPair();
		}

		@Override
		public void close() throws IOException
		{
			lineReader.close();
		}

		@Override
		public LongWritable getCurrentKey() 
		{
			return lineReader.getCurrentKey();
		}

		@Override
		public ReadPair getCurrentValue() 
		{
			return value;
		}

		@Override
		public float getProgress() throws IOException 
		{
			return lineReader.getProgress();
		}

		@Override
		public boolean nextKeyValue() throws IOException
		{
			if (lineReader.nextKeyValue())
			{
				Text line = lineReader.getCurrentValue();
				value.clear();

				try
				{
					TextSamMapping mapping = new TextSamMapping(line);
					if (mapping.isRead2())
						value.setRead2(mapping);
					else // anything that's not explicitly labelled as "read 2" goes in as read 1.
						value.setRead1(mapping);

					return true;
				}
				catch (FormatException e) {
					throw new FormatException(e.getMessage() + ". File: " + split.getPath() + "; Position: " + lineReader.getCurrentKey().get());
				}
			}
			else
				return false;
		}
	}

	@Override
	public RecordReader<LongWritable,ReadPair>	createRecordReader(InputSplit split, TaskAttemptContext context)
	{
		FileSplit fsplit = (FileSplit)split;

		if (fsplit.getStart() > 0 && !isSplitable(context, fsplit.getPath()))
			throw new RuntimeException("Trying to split non-splittable file " + fsplit.getPath());

		return new SamRecordReader();
	}

	@Override
	protected boolean isSplitable(JobContext context, Path file)
	{
		CompressionCodec codec = new CompressionCodecFactory(context.getConfiguration()).getCodec(file);
		return codec == null;
	}
}
