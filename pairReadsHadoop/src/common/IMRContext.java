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

/**
 * An interface that covers the basic functionality we require from the Map/Reduce Context objects.
 * By programming to this interface instead of the bare Hadoop Mapper.Context and Reducer.Context
 * objects we should be able to more easily unit-test the Mappers and Reducers.
 */
public interface IMRContext<KEYOUT, VALUEOUT>
{
	public void progress();
	public void setStatus(String msg);
	public void write(KEYOUT key, VALUEOUT value) throws java.io.IOException, InterruptedException;

	public void increment(Enum<?> counterName, long value);
	public void increment(String groupName, String counterName, long value);
}
