//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>
#include <algorithm>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#ifdef _DEBUG_LOG
#define DEBUG_PRINT(x) do {std::cerr << timestamp() << " " << x << std::endl;} \
					   while(0)
#else
#define DEBUG_PRINT(x)
#endif

inline std::string timestamp(const char* format = "[%H:%M:%S]")
{
	std::time_t t = std::time(0);
	char cstr[128];
	std::strftime(cstr, sizeof(cstr), format, std::localtime(&t));
	return cstr;
}

//taken from http://stackoverflow.com/questions/14826893
template <typename FwdIt> class adjacent_iterator 
{
public:
    adjacent_iterator(FwdIt first, FwdIt last)
        : m_first(first), m_next(first == last ? first : std::next(first)) { }

    bool operator!=(const adjacent_iterator& other) const 
	{
        return m_next != other.m_next; // NOT m_first!
    }

    adjacent_iterator& operator++() 
	{
        ++m_first;
        ++m_next;
        return *this;
    }

    typedef typename std::iterator_traits<FwdIt>::reference Ref;
    typedef std::pair<Ref, Ref> Pair;

    Pair operator*() const 
	{
        return Pair(*m_first, *m_next); // NOT std::make_pair()!
    }

private:
    FwdIt m_first;
    FwdIt m_next;
};

template <typename FwdIt> class adjacent_range 
{
public:
    adjacent_range(FwdIt first, FwdIt last)
        : m_first(first), m_last(last) { }

    adjacent_iterator<FwdIt> begin() const 
	{
        return adjacent_iterator<FwdIt>(m_first, m_last);
    }

    adjacent_iterator<FwdIt> end() const 
	{
        return adjacent_iterator<FwdIt>(m_last, m_last);
    }

private:
    FwdIt m_first;
    FwdIt m_last;
};

template <typename C> auto make_adjacent_range(C& c) -> 
	adjacent_range<decltype(c.begin())> 
{
    return adjacent_range<decltype(c.begin())>(c.begin(), c.end());
}

template <class C, class V>
bool contains(C container, V value)
{
	return std::count(container.begin(), container.end(), value);
}

template <class T>
void vecRemove(std::vector<T>& vec, const T& val)
{
	vec.erase(std::remove(vec.begin(), vec.end(), val), std::end(vec));
}

inline bool makeDirectory(const std::string& name)
{
#ifdef _WIN32
	int ret = _mkdir(name.c_str());
#else
	int ret = mkdir(name.c_str(), 0777);
#endif
	return !ret;
}
