#ifndef DIGITALSNOW_DEFORMATIONS_GRAYSCALEMAPCAST_HPP
#define DIGITALSNOW_DEFORMATIONS_GRAYSCALEMAPCAST_HPP

// Grayscale map returning value casted to specified type
template <typename T, typename TCast = unsigned char>
class GrayscaleMapCast
{
public:
    GrayscaleMapCast() : m_min(0), m_max(255) {}
    GrayscaleMapCast(T min_value, T max_value) : m_min(min_value), m_max(max_value) {}
    
    inline
    unsigned char operator() (T value) const
    {
        return static_cast<TCast>( (255*(value - m_min)) / (m_max - m_min) );
    }

private:
    T m_min, m_max;
};

#endif // DIGITALSNOW_DEFORMATIONS_GRAYSCALEMAPCAST_HPP
