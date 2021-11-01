1. 实现单帧降噪,新增求单个像素降噪结果函数FilterPixel
2. 实现两帧间的投影
3. 实现两帧间的累积，新增单个像素Clamp函数ClampAcc

ps:当场景为pink-room时，m_sigmaColor改为10