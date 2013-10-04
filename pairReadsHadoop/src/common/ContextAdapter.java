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

import common.IMRContext;

import org.apache.hadoop.mapreduce.Counter;

public class ContextAdapter<K, V> implements IMRContext<K, V>
{
	private org.apache.hadoop.mapreduce.TaskInputOutputContext<?, ?, K, V> hadoopContext;

	public ContextAdapter(org.apache.hadoop.mapreduce.TaskInputOutputContext<?, ?, K, V> m)
	{
		hadoopContext = m;
	}

	public void progress()
	{
		hadoopContext.progress();
	}

	public void setStatus(String msg)
	{
		hadoopContext.setStatus(msg);
	}

	public void write(K k, V v) throws java.io.IOException, InterruptedException
	{
		hadoopContext.write(k, v);
	}

	public void increment(Enum<?> counterName, long value)
	{
		hadoopContext.getCounter(counterName).increment(value);
	}

	public void increment(String groupName, String counterName, long value)
	{
		hadoopContext.getCounter(groupName, counterName).increment(value);
	}
}
