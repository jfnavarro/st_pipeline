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

import java.io.DataOutputStream;
import java.io.IOException;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.util.ReflectionUtils;

public class OutputStreamFactory
{
	protected TaskAttemptContext context;

	public OutputStreamFactory(TaskAttemptContext context)
	{
		this.context = context;
	}

	public DataOutputStream makeStream(Path path) throws IOException
	{
		Configuration conf = context.getConfiguration();
		boolean isCompressed = FileOutputFormat.getCompressOutput(context);

		CompressionCodec codec = null;
		String extension = "";

		if (isCompressed)
		{
			Class<? extends CompressionCodec> codecClass = FileOutputFormat.getOutputCompressorClass(context, GzipCodec.class);
			codec = (CompressionCodec) ReflectionUtils.newInstance(codecClass, conf);
			extension = codec.getDefaultExtension();
		}

		FileSystem fs = path.getFileSystem(conf);

		DataOutputStream output;

		if (isCompressed)
		{
			FSDataOutputStream fileOut = fs.create(path, false);
			output = new DataOutputStream(codec.createOutputStream(fileOut));
		}
		else
			output = fs.create(path, false);

		return output;
	}
}
