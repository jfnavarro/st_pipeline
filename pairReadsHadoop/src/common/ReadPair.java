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

import java.util.Iterator;
import java.util.NoSuchElementException;

public class ReadPair implements Iterable<AbstractTaggedMapping>
{
	protected AbstractTaggedMapping read1;
	protected AbstractTaggedMapping read2;

	public ReadPair()
	{
		read1 = read2 = null;
	}

	public ReadPair(AbstractTaggedMapping r1, AbstractTaggedMapping r2)
	{
		/*
		if (r1 != null && r1.isPaired() && !r1.isRead1())
			throw new IllegalArgumentException("first read in pair is marked as paired but not as read 1");
		if (r2 != null && r2.isPaired() && !r2.isRead2())
			throw new IllegalArgumentException("second read in pair is marked as paired but not as read 2");
			*/
		read1 = r1;
		read2 = r2;
	}

	public void clear()
	{
		read1 = read2 = null;
	}

	public AbstractTaggedMapping getAnyRead()
	{
		if (read1 != null)
			return read1;
		else
			return read2;
	}

	public String getName()
	{
		AbstractTaggedMapping read = getAnyRead();
		if (read == null)
			return null;

		String name = read.getName();
		// trim the trailing /n, if any
		if (name.length() >= 2)
		{
			int last = name.length() - 1;
			if (name.charAt(last - 1) == '/' && name.charAt(last) >= '0' && name.charAt(last) <= '9')
				name = name.substring(0, name.length() - 2);
		}

		return name;
	}

	public AbstractTaggedMapping getRead1() { return read1; }
	public AbstractTaggedMapping getRead2() { return read2; }

	public void setRead1(AbstractTaggedMapping map) { read1 = map; }
	public void setRead2(AbstractTaggedMapping map) { read2 = map; }


	///////////////// iterator ///////////////////

	protected static class ReadIt implements Iterator<AbstractTaggedMapping>
	{
		/*
		 * class invariant:
		 * state always indicates the next non-null read to return.
		 * Once returned, the read variables are set to null, then
		 * the iterator is advanced.
		 */
		private AbstractTaggedMapping r1, r2;
		private enum IteratorState { One, Two, End }; // next read to check
		private IteratorState state;

		public ReadIt(AbstractTaggedMapping a, AbstractTaggedMapping b)
		{
			r1 = a;
			r2 = b;
			advance();
		}

		private void advance()
		{
			if (r1 == null)
			{
				if (r2 == null)
					state = IteratorState.End;
				else
					state = IteratorState.Two;
			}
			else
				state = IteratorState.One;
		}


		public boolean hasNext()
		{
			return state != IteratorState.End;
		}

		public AbstractTaggedMapping next()
		{
			AbstractTaggedMapping retval = null;
			if (state == IteratorState.One)
			{
				retval = r1;
				r1 = null;
			}
			else if (state == IteratorState.Two)
			{
				retval = r2;
				r2 = null;
			}
			else
				throw new NoSuchElementException();

			advance();
			return retval;
		}

		public void remove()
		{
			throw new UnsupportedOperationException();
		}
	}

	public Iterator<AbstractTaggedMapping> iterator()
	{
		return new ReadIt(read1, read2);
	}
}
