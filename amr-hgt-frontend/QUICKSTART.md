# amr-hgt-frontend

Pre-built Vite + React project — no `npm create vite` needed.

## Setup (run once)

```bash
cd amr-hgt-frontend
npm install
```

## Run

```bash
npm run dev
```

Open the URL it prints — usually http://localhost:5173

## Before running

Make sure the backend is also running, in a separate terminal:

```bash
cd ../amr-hgt-tool/backend
python app.py
```

It should print `Running on http://127.0.0.1:5050`. Leave both terminals
open at the same time — this frontend (5173) is the page you open in your
browser; it calls the backend (5050) when you click "Analyze".

## If you update the dashboard code later

Just replace `src/App.jsx` with the new version — no need to redo `npm install`
unless dependencies changed.
