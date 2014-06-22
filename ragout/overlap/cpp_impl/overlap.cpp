//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>

#include "overlap.h"
#include "fasta.h"
#include "suffix_array.h"

namespace {

struct coords
{
	long long numOfString;
	long long numOfSuffix;
	coords(long long nOfString, long long nOfSuffix):
		numOfString(nOfString), numOfSuffix(nOfSuffix)
	{}

	coords():
		numOfString(0), numOfSuffix(0)
	{}

};

struct FMIndex
{
	static const int ALPHABET = 255;
	std::unordered_map<char, std::vector<long long>> occ_;
	std::unordered_map<char, long long> c_;
	std::vector<char> BWT_;
	std::vector<coords> gsa_;
	std::vector<int> lexOrder_;
	FMIndex()
	{}

	FMIndex(std::vector<char>& BWT, std::vector<coords>& gsa) :
		BWT_(BWT), gsa_(gsa)
	{
		for (size_t i = 0; i < gsa.size(); ++i)
		{
			if (gsa[i].numOfSuffix == 0)
			{
				lexOrder_.push_back(gsa[i].numOfString);
			}
		}
		std::vector<long long> temp(ALPHABET);
		for (size_t i = 0; i < BWT_.size(); ++i)
		{
			temp[BWT_[i] + 1]++;
		}

		for (size_t i = 0; i < temp.size(); ++i)
		{
			if (temp[i] != 0)
			{
				occ_[(char)i - 1].resize(BWT_.size());
			}
		}
		for (size_t i = 1; i < temp.size(); ++i)
		{
			temp[i] += temp[i - 1];
			c_[(char)i] = temp[i];

		} 
		occ_[BWT_[0]][0] = 1;
		for (auto it : occ_)
		{
			for (size_t i = 1; i < BWT_.size(); ++i)
			{
				if (it.first == BWT_[i]) 
				{
					occ_[it.first][i]++;
				}
				occ_[it.first][i] += occ_[it.first][i - 1];
			}
		}
	}

	std::vector<Overlap> findOverlaps(std::vector<FastaRecord>& fRecords,
									  int minOverlap)
	{
		std::vector<Overlap> overlaps;
		std::cerr << "\tOverapping\n";
		auto cmp = [] (Overlap* o1, Overlap* o2) {return o1->size < o2->size;};
		for (size_t i = 0; i < fRecords.size(); ++i)
		{
			std::vector<Overlap> contigOverlaps;
			std::unordered_map<FastaRecord*, std::vector<Overlap*>> ovlpTable;

			findOccures(i, minOverlap, fRecords, contigOverlaps);
			for (Overlap& ovlp : contigOverlaps)
			{
				ovlpTable[ovlp.nextContig].push_back(&ovlp);
			}
			for (auto tablePair : ovlpTable)
			{
				overlaps.push_back(**std::max_element(tablePair.second.begin(),
													  tablePair.second.end(),
													  cmp));
			}
		}
		return overlaps;
	}

	void findOccures(int indexInRecord, int minOverlap,
					 std::vector<FastaRecord>& fRecords,
					 std::vector<Overlap>& overlapsOut)
	{
		std::string read = fRecords[indexInRecord].sequence_;
		if (read[read.length() - 1] == 'N') return;
		long long l = c_[read[read.length() - 1]];
		long long r = c_[read[read.length() - 1] + 1] - 1;
		int i = read.length() - 2;
		
		while(l + 1 <= r && i >= 0)
		{
			if (read[i] == 'N') return;
			l = c_[read[i]] + occ_[read[i]][l - 1];
			r = c_[read[i]] + occ_[read[i]][r] - 1;

			if ((int)read.length() - i >= minOverlap)
			{
				long long  newl = c_['$'] + occ_['$'][l - 1];
				long long newr = c_['$'] + occ_['$'][r] - 1;
				if (newl <= newr)
				{
					for (int j = newl; j <= newr; ++j)
					{
						if (lexOrder_[j] != indexInRecord)
						{
							overlapsOut.push_back(Overlap(&fRecords[indexInRecord],
														  &fRecords[lexOrder_[j]],
														  read.length() - i));
						}
					}
				}
			}
			--i;
		}
	}

};

//transform suffix array at form (i,j) -> i - NumberOfString, 
//j - numberOfSuffix of that string
void getGSA(std::vector<long long>& SA, std::string& superText,
			std::vector<coords>& GSA)
{

	std::vector<coords> preGSA;
	preGSA.reserve(SA.size());
	GSA.resize(SA.size());
	long long numberof$ = 0;
	long long numSuff = 0;
	for (size_t i = 0; i < superText.length(); ++i)
	{
		preGSA.push_back(coords(numberof$, numSuff));
		numSuff++;
		if (superText[i] == '$')
		{
			numberof$++;
			numSuff = 0;
		}
	}
	for (size_t i = 0; i < superText.length(); ++i)
	{
		GSA[i] = preGSA[SA[i]];
	}
}

void getBwt(std::vector<coords>& GSA, std::vector<FastaRecord>& fRecords,
			std::vector<char>& BWT)
{
	BWT.resize(GSA.size());
	for (size_t i = 0; i < GSA.size(); ++i)
	{
		if (GSA[i].numOfSuffix != 0)
		{
			BWT[i] = fRecords[GSA[i].numOfString].sequence_[GSA[i].numOfSuffix - 1];
			
			assert(std::string("ACGTN").find(BWT[i]) != std::string::npos);
		}
		else
		{
			BWT[i] = '$';
		}
	}
}

//concatenate everything with $
std::string getSuperString(std::vector<FastaRecord>& fRecords,
						   int maxOverlapLength)
{
	std::string superString = "";
	long long length = 0;
	for (size_t i = 0; i < fRecords.size(); ++i)
	{
		if ((int)fRecords[i].sequence_.length() <= 2 * maxOverlapLength)
		{
			length += fRecords[i].sequence_.length() + 1;
		}
		else
		{
			length += 2 * maxOverlapLength + 1;
		}
	}
	superString.reserve(length);
	//int eliminated = 0;
	for (size_t i = 0; i < fRecords.size(); ++i)
	{
		if ((int)fRecords[i].sequence_.length() <= 2 * maxOverlapLength)
		{
			superString.append(fRecords[i].sequence_);
			superString.push_back('$');
		}
		else
		{
			superString += fRecords[i].sequence_.substr(0, maxOverlapLength);
			superString += fRecords[i].sequence_.substr(fRecords[i].sequence_.length() -
														maxOverlapLength, 
														maxOverlapLength);
			superString += "$";
			fRecords[i].sequence_ = fRecords[i].sequence_.substr(0, maxOverlapLength) + 
									fRecords[i].sequence_.substr(fRecords[i].sequence_.length() -
																 maxOverlapLength, maxOverlapLength);
		}
	}
	return superString;
}

FMIndex getFM(std::vector<FastaRecord>& fRecords, int maxOverlapLength)
{
	std::cerr << "\tBuilding FM-index\n";
	std::string superString = getSuperString(fRecords, maxOverlapLength);
	unsigned char* text = new unsigned char[superString.size()];
	long long* SA = new long long[superString.size()];
	std::copy(superString.begin(), superString.end(), text);

	sais(text, SA, superString.length());
	std::vector<long long> vSA(superString.length());
	std::copy(SA, SA + superString.length(), vSA.begin());
	delete[] SA;
	delete[] text;

	std::vector<coords> gsa;
	getGSA(vSA, superString, gsa);
	superString.clear();
	vSA.clear();

	std::vector<char> BWT;
	getBwt(gsa, fRecords, BWT);
	FMIndex FM(BWT, gsa);
	return FM;
}

} //end anonymous namespace

std::vector<Overlap> getOverlaps(std::vector<FastaRecord>& contigs, 
								 int minOverlap, int maxOverlap)
{
	FMIndex FM = getFM(contigs, maxOverlap);
	return FM.findOverlaps(contigs, minOverlap);
}
