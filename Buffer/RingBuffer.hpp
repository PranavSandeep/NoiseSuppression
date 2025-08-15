#include <cassert>
#include <atomic>
#include <iostream>
#include <vector>

template <class T>
class RingBuffer
{
    public:
    /**
     *
     * @param capacity The total size of the ring buffer. Must be a power of 2.
     */
    explicit RingBuffer(int capacity)
        : capacity(capacity), buffer(new T[capacity])
    {
        assert((capacity & (capacity-1)) == 0); //Checking if power of 2.
        mask = capacity - 1;
        head_.store(0, std::memory_order_relaxed);
        tail_.store(0, std::memory_order_relaxed);
    }

    void addFromVector(std::vector<T> const &vector)
    {
        size_t size = vector.size();
        size_t tail = tail_.load(std::memory_order_relaxed);
        size_t head = head_.load(std::memory_order_acquire);

        if (size > capacity)
            size = capacity;

        size_t start = tail & (mask);
        size_t space_at_end = capacity - start;
        if (space_at_end > size)
            std::memcpy(buffer.get() + start, vector.data(), size * sizeof(T));
        else
        {
            size_t first = space_at_end;
            size_t second = size-first;
            std::memcpy(buffer.get() + start, vector.data(), first * sizeof(T));
            std::memcpy(buffer.get(), vector.data()+first, second * sizeof(T));
        }

        tail_.store(tail + size, std::memory_order_release);

    }

    T readFromBuffer() const
    {
        return buffer[head_.load(std::memory_order_relaxed)];
    }

    void seek(int offset)
    {
        head_.store((head_+offset) & mask, std::memory_order_relaxed);
    }

    [[nodiscard]] int getHead() const
    {
        return head_.load(std::memory_order_relaxed);
    }


    ~RingBuffer() = default;

    private:
    size_t mask;
    size_t capacity;
    std::unique_ptr<T[]> buffer;
    std::atomic<size_t> head_;
    std::atomic<size_t> tail_;
};
