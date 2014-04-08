//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include <stdexcept>

namespace
{
	class ParseException : public std::runtime_error 
	{
	public:
		ParseException(const std::string & what):
			std::runtime_error(what)
		{}
	};

	char complementSymbol(char c)
	{
		switch (c)
		{
			case 'A': return 'T';
			case 'C': return 'G';
			case 'G': return 'C';
			case 'T': return 'A';
			case 'a': return 't';
			case 'c': return 'g';
			case 'g': return 'c';
			case 't': return 'a';
			default: return c;
		}
	}

	void reverseComplement(const std::string& dnaIn, std::string& dnaOut)
	{
		dnaOut.clear();
		for (size_t i = dnaIn.length(); i > 0; --i)
		{
			dnaOut += complementSymbol(dnaIn[i - 1]);
		}
	}
}

size_t FastaReader::GetSequences(std::vector<FastaRecord> & record)
{
	std::string buffer;
	std::string sequence;
	std::string header;
	int line = 1;
	size_t seqId = record.size();

	try
	{
		while(!inputStream_.eof())
		{
			std::getline(inputStream_, buffer, '\n');
			if (*buffer.rbegin() == '\r') buffer.erase(buffer.size() - 1);
			if (buffer.empty()) continue;

			if (buffer[0] == '>')
			{
				if (!header.empty())
				{
					if (sequence.empty()) throw ParseException("empty sequence");

					record.push_back(FastaRecord(sequence, header, seqId));
					++seqId;
					sequence.clear();
					header.clear();
				}
				ValidateHeader(buffer);
				header = buffer;
			}
			else
			{
				ValidateSequence(buffer);
				sequence += buffer;
			}

			++line;
		}
		
		if (sequence.empty()) throw ParseException("empty sequence");
		record.push_back(FastaRecord(sequence, header, seqId));
	}
	catch (ParseException & e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName_ << " on line " << line << ": " << e.what();
		throw std::runtime_error(ss.str());
	}

	return record.size();
}

size_t FastaReader::GetSequencesWithComplements(std::vector<FastaRecord>& records)
{
	this->GetSequences(records);
	std::vector<FastaRecord> complements;
	for (auto &record : records)
	{
		std::string header = "-" + record.description_;
		record.description_ = "+" + record.description_;
		std::string revComplement;
		reverseComplement(record.sequence_, revComplement);
		complements.push_back(FastaRecord(revComplement, header, record.id_));
	}
	std::copy(complements.begin(), complements.end(), std::back_inserter(records));
	return records.size();
}

void FastaReader::ValidateHeader(std::string & header)
{
	size_t delim = header.find(' ');
	if (delim == std::string::npos)
	{
		delim = header.length() - 1;
	}
	else
	{
		--delim;
	}

	header = header.substr(1, delim);
	if (header.empty()) throw ParseException("empty header");
}

void FastaReader::ValidateSequence(std::string & sequence)
{
	const std::string VALID_CHARS = "ACGTURYKMSWBDHWNX-";
	for (size_t i = 0; i < sequence.length(); ++i)
	{
		char orig = sequence[i];
		sequence[i] = toupper(sequence[i]);
		if (VALID_CHARS.find(sequence[i]) == std::string::npos) 
		{
			throw ParseException((std::string("illegal character: ") + orig).c_str());
		}
	}
}

bool FastaReader::IsOk() const
{
	return inputStream_.good();
}
