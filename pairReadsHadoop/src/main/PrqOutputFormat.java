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

import common.AbstractTaggedMapping;
import common.OutputStreamFactory;
import common.ReadPair;

import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

/**
 * Output format for the prq format.
 *
 */
public class PrqOutputFormat extends FileOutputFormat<Text, ReadPair>
{
	public static class PrqRecordWriter extends RecordWriter<Text,ReadPair>
	{
		private DataOutputStream os;

		public PrqRecordWriter(DataOutputStream stream)
		{
			os = stream;
		}

		/**
		 * Writes one mapping.
		 *
		 * Writes \t sequence \t quality.  If map is null it just writes
		 * the delimiters with empty fields.
		 * 
		 * Not that it doesn't write a trailing delimiter.
		 */
		private void writeMapping(AbstractTaggedMapping map) throws IOException
		{
			os.writeByte('\t');
			if (map != null)
			{
				ByteBuffer buffer = map.getSequence();
				os.write(buffer.array(), buffer.position(), map.getLength());
				os.writeByte('\t');
				buffer = map.getBaseQualities();
				os.write(buffer.array(), buffer.position(), map.getLength());
			}
			else
				os.writeByte('\t');
		}

		public void write(Text key, ReadPair pair) throws IOException
		{
			if (key != null)
				os.write(key.getBytes(), 0, key.getLength());
			else
				os.writeBytes(pair.getName());

			writeMapping(pair.getRead1());
			writeMapping(pair.getRead2());
			os.writeByte('\n');
		}

		public void close(TaskAttemptContext context) throws IOException
		{
			os.close();
		}
	}

  public RecordWriter<Text,ReadPair> getRecordWriter(TaskAttemptContext task)
	  throws IOException
	{
		DataOutputStream os = new OutputStreamFactory(task).makeStream(this.getDefaultWorkFile(task, ""));
		return new PrqRecordWriter(os);
	}
}
