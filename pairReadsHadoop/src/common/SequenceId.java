/*
 * Copyright (C) 2011-2012 CRS4.
 *
 * This file is part of Seal.
 *
 * Seal is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Seal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Seal.  If not, see <http://www.gnu.org/licenses/>.
 */

package common;

import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.WritableComparable;
import org.apache.hadoop.io.WritableComparator;
import org.apache.hadoop.io.WritableUtils;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class SequenceId implements WritableComparable<SequenceId>
{
	private String location;
	private byte read;

	public SequenceId()
	{
		location = "";
		read = 127;
	}

	public SequenceId(String location, int read)
	{
		set(location, read);
	}

	public void set(String location, int read)
	{
		if (read <= 0 || read > 127)
			throw new IllegalArgumentException("invalid read number " + read);

		this.location = location;
		this.read = (byte)read;
	}

	public String getLocation() { return location; }
	public int getRead() { return read; }

	/**
	 * Read the two fields.
	 * Encoded as:
	 * 	location length: Vint
	 * 	location: UTF bytes
	 * 	read: byte
	 */
	@Override
	public void readFields(DataInput in) throws IOException
	{
		location = Text.readString(in);
		read = in.readByte();
	}

	@Override
	public void write(DataOutput out) throws IOException
	{
		Text.writeString(out, location);
		out.writeByte(read);
	}

	@Override
	public int hashCode()
	{
		return location.hashCode() + 157*read;
	}

	@Override
	public boolean equals(Object right)
	{
		if (right instanceof SequenceId)
		{
			SequenceId r = (SequenceId) right;
			return r.location.equals(location) && r.read == read;
		}
		else
			return false;
	}

	@Override
	public int compareTo(SequenceId o)
	{
		int locationComparison = location.compareTo(o.location);
		if (locationComparison != 0)
			return locationComparison;
		else if (read != o.read)
			return read < o.read ? -1 : 1;
		else
			return 0;
	}

	@Override
	public String toString()
	{
		return "(" + location + "," + read + ")";
	}

	/** A Comparator that compares serialized SequenceId. */
	public static class Comparator extends WritableComparator
	{
		public Comparator()
		{
			super(SequenceId.class);
		}

		public int compare(byte[] b1, int s1, int l1, byte[] b2, int s2, int l2)
		{
			int sizeVint1 = WritableUtils.decodeVIntSize(b1[s1]);
			int sizeVint2 = WritableUtils.decodeVIntSize(b2[s2]);
			return compareBytes(b1, s1+sizeVint1, l1-sizeVint1, b2, s2+sizeVint2, l2-sizeVint2);
		}
	}

	static {   // register this comparator
		WritableComparator.define(SequenceId.class, new Comparator());
	}
}
