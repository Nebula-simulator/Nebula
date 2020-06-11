#ifndef __OUTPUT_STREAM_H_
#define __OUTPUT_STREAM_H_

#include <iostream>
#include <fstream>
#include <mutex>
#include <algorithm>
#include <array>
#include <assert.h>

/**
 * \brief Thread-safe binary output stream.
 *
 * Will write either to files on disk or to the standard output stream.
 */

class output_stream
{
public:
	/**
	 * \brief Constructor.
	 *
	 * If the special file name `stdout` is provided, this class will write to
	 * the standard output stream instead of a file.
	 *
	 * \param filename File name to write to, or `stdout` for the standard
	 *                 output stream.
	 */
	explicit output_stream(std::string const & filename)
	{
		if (filename == "stdout")
		{
			file = &std::cout;
			own_file = false;
		}
		else
		{
			file = new std::ofstream(filename, std::ofstream::binary);
			own_file = true;
		}
	}

	~output_stream()
	{
		file->flush();
		if (own_file)
			delete file;
	}

	/// Write data in a buffer to the file.
	void write(char const * data, size_t size)
	{
		std::lock_guard<std::mutex> lock(mutex);
		file->write(data, size);
	}

private:
	std::mutex mutex;
	std::ostream* file = nullptr;
	bool own_file = false;

	output_stream(output_stream const &) = delete;
	output_stream& operator=(output_stream const &) = delete;
	output_stream(output_stream&&) = delete;
	output_stream& operator=(output_stream&&) = delete;
};

/**
 * Buffer to be used with ::output_stream.
 *
 * This class is not, in itself, thread-safe. It is intended for thread-local
 * use. This class will automatically flush its data to a thread-safe
 * ::output_stream when the buffer is full.
 *
 * If data with multiple data types need to be written and the order needs
 * preserving, there are two options.
 *   1. Ensure that the buffer size is exactly a multiple of each "data row",
 *      e.g. if you need to write 3 `float`s and 1 `int`, set the buffer size to
 *      any multiple of `(3*sizeof(float) + 1*sizeof(int))`
 *   2. Merge the contents in a temporary `char` buffer array yourself.
 */
class output_buffer
{
public:
	/**
	 * \brief Constructor.
	 *
	 * \param stream ::output_stream this buffer prints to
	 * \param size   Number of bytes in the buffer
	 */
	output_buffer(output_stream& stream, size_t size)
		: buffer(size), position(0), outstream(stream)
	{}

	/**
	 * \brief Add data to the buffer.
	 *
	 * The buffer must be large enough to hold the data. If there is not enough
	 * space in the buffer, the data is flushed to the output stream.
	 */
	void add(char const * data, size_t size)
	{
		assert(size < buffer.size());
		if (position + size > buffer.size())
			flush();
		std::copy_n(data, size, buffer.data() + position);
		position += size;
	}

	/**
	 * \brief Add data to the buffer.
	 *
	 * Data is passed as a `std::array`. All data in the array is put into the
	 * buffer. The buffer must be large enough to hold the array data at once.
	 * If there is not enough space in the buffer, the data is flushed to the
	 * output stream.
	 */
	template<typename T, size_t N>
	void add(std::array<T, N> const & data)
	{
		add(reinterpret_cast<char const *>(data.data()), N*sizeof(T));
	}

	/// Send the buffer contents to the output stream, and clear it.
	void flush()
	{
		outstream.write(buffer.data(), position);
		position = 0;
	}

private:
	std::vector<char> buffer;
	size_t position;
	output_stream& outstream;
};

#endif // __OUTPUT_STREAM_H_
