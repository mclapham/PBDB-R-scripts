import express from 'express';
import cors from 'cors';
import fetch from 'node-fetch';

const app = express();
app.use(cors());
app.use(express.json());
app.use(express.urlencoded({ extended: true }));

app.post('/pbdb-update', async (req, res) => {
  const { session_id, ...payload } = req.body;

  if (!session_id) {
    return res.status(400).json({ error: 'Missing session_id parameter' });
  }

  const params = new URLSearchParams(payload);

  try {
    const pbdbRes = await fetch('https://paleobiodb.org/data1.2/colls/update.json', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Cookie': `session_id=${session_id}`
      },
      body: params.toString()
    });

    const data = await pbdbRes.json().catch(() => null);
    res.status(pbdbRes.status).json(data);
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => console.log(`Proxy listening on port ${PORT}`));
