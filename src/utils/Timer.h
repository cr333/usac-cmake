#ifndef TIMER_H
#define TIMER_H

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/time.h>
#endif

namespace Timer
{
	inline double getTimestampInSeconds()
	{
#if _WIN32

		LARGE_INTEGER li_time, li_freq;
		QueryPerformanceCounter(  &li_time);
		QueryPerformanceFrequency(&li_freq);

		return double(li_time.QuadPart) / double(li_freq.QuadPart);

#else

		struct timeval timestamp;
		gettimeofday(&timestamp, NULL);

		return timestamp.tv_sec + timestamp.tv_usec / 1000000.0;

#endif
	}
}

#endif
