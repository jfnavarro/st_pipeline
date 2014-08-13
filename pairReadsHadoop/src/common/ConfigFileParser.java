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

import java.io.IOException;
import java.io.File;
import java.io.LineNumberReader;
import java.io.Reader;

import java.util.regex.*;
import java.util.Properties;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Collection;


public class ConfigFileParser
{
	// Line patterns.  Lines are trimmed of spaces before being matched.
	// Section name
	private static final Pattern PSectionTitle = Pattern.compile("^\\[ *(\\w+) *\\] *$");
	// key-value.  Key is delimited by : or =.  Spaces before the value are removed.  Since
	// the line is trimmed, spaces at the end of the value will also be removed.
	// Except for this latest point, this should match the behaviour of the Python
	// ConfigParser module.
	private static final Pattern PKvLine = Pattern.compile("^ *([^:=]+)[:=] *(.*) *$");

	// Settings in the default section are inherited by all other sections.
	// It must be named DEFAULT for compatibility.
	private static final String DefaultSectionName = "DEFAULT";

	private Properties defaultProperties;
	private HashMap<String, Properties> sections;

	private enum TokenType {
		SectionName,
		KeyValue
	}

	public static class KvPair
	{
		private String key;
		private String value;

		public KvPair(String k, String v) {
			key = k;
			value = v;
		}

		public String getKey() { return key; }
		public String getValue() { return value; }

		public boolean equals(Object other) {
			if (other instanceof KvPair)
			{
				KvPair otherPair = (KvPair)other;
				if (key == null && value == null)
					return otherPair.key == null && otherPair.value == null;
				else if (key == null && value != null)
					return (key == otherPair.key) && value.equals(otherPair.value);
				else if (value == null && key != null)
					return (value == otherPair.value) && key.equals(otherPair.key);
				else
					return key.equals(otherPair.key) && value.equals(otherPair.value);
			}
			else
				return false;
		}

		public int hashCode() {
			int k = (key == null) ? 0 : key.hashCode();
			int v = (value == null) ? 0 : value.hashCode();
			return k ^ v;
		}

		public String toString() {
		 	return ((key == null) ? "''" : key) + " => " + ((value == null) ? "''" : value);
		}
	}

	public static class KvIterator implements Iterator<KvPair>
	{
		private Iterator<String> namesIt;
		private Properties properties;

		public KvIterator(Properties p)
		{
			properties = p;
			namesIt = p.stringPropertyNames().iterator();
		}

		public boolean hasNext()
		{
			return namesIt.hasNext();
		}

		public KvPair next()
		{
			String propName = namesIt.next();
			return new KvPair(propName, properties.getProperty(propName));
		}

		public void remove()
		{
			throw new UnsupportedOperationException();
		}
	}

	// in here we maintain the current scanner state and token
	private LineNumberReader input;
	private TokenType currentToken;

	private String currentSectionName;
	private KvPair tokenKV;

	public ConfigFileParser()
	{
		defaultProperties = new Properties();
		sections = new HashMap<String, Properties>();
		currentSectionName = null;;
	}

	/**
	 * Reads the next high-level token from the input file.
	 * Advances the reader to the next high-level token, either a section title or
	 * a key-value assignment.
	 *
	 * @return true if the next token was read; false if the token was not advanced because there's no more to read.
	 * @exception IOException Something bad happened while reading the file.
	 * @exception FormatException An invalid line was read from the config file.
	 */
	private boolean advanceToken() throws IOException, FormatException
	{
		String line = null;
		boolean done = false;
		do
		{
			line = input.readLine();
			if (line != null) // readLine returns null on EOF
			{
				line = line.trim();
				if (line.isEmpty() || line.charAt(0) == '#' || line.charAt(0) == ';')
				{
					// comment line.  Skip
				}
				else
				{
					Matcher titleMatcher = PSectionTitle.matcher(line);
					Matcher kvMatcher = PKvLine.matcher(line);
					if (titleMatcher.matches())
					{
						currentToken = TokenType.SectionName;
						currentSectionName = titleMatcher.group(1);
						done = true;
					}
					else if (kvMatcher.matches())
					{
						currentToken = TokenType.KeyValue;
						tokenKV = new KvPair(kvMatcher.group(1).trim(), kvMatcher.group(2).trim());
						done = true;
					}
					else
						throw new FormatException("Invalid config line (" + input.getLineNumber() + "): " + line);
				}
			}
			else
				done = true; // EOF
		}
		while (!done);

		return line != null;
	}

	public void load(Reader cfgFile) throws IOException, FormatException
	{
		input = new LineNumberReader(cfgFile);

		// re-initialize data structures to support reloading
		defaultProperties = new Properties();
		sections = new HashMap<String, Properties>();
		currentSectionName = null;

		while (advanceToken())
		{
			if (TokenType.KeyValue == currentToken)
			{
				Properties p;
				if (currentSectionName == null || currentSectionName.equalsIgnoreCase(DefaultSectionName))
					p = defaultProperties;
				else
					p = sections.get(currentSectionName);

				p.put(tokenKV.getKey(), tokenKV.getValue());
			}
			else if (TokenType.SectionName == currentToken)
			{
				sections.put(currentSectionName, new Properties(defaultProperties));
			}
			else
				throw new RuntimeException("Unknown token type " + currentToken + ". Please file a bug report");
		}
	}

	public Collection<String> getSectionNames()
	{
		return sections.keySet();
	}

	public boolean hasSection(String sectionName)
	{
		if (sectionName != null && sectionName.equalsIgnoreCase(DefaultSectionName))
			return true;
		else
			return sections.containsKey(sectionName);
	}

	public Iterator<KvPair> getSectionIterator(String sectionName)
	{
		Properties p = sections.get(sectionName);
		if (p == null)
			p = new Properties(defaultProperties);

		return new KvIterator(p);
	}

	public String getValue(String sectionName, String key)
	{
		Properties section = sections.get(sectionName);
		if (section != null)
			return section.getProperty(key);
		else
			return defaultProperties.getProperty(key);
	}
}
